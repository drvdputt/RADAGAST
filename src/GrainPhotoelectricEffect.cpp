#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "Testing.h"

#include <cmath>
#include <iomanip>
#include "GrainPhotoelectricEffect.h"

//#define EXACTGRID

using namespace std;

namespace
{
vector<double> _filelambdav, _fileav;
vector<vector<double>> _Qabsvv, _Qscavv, _asymmparvv;
}

void GrainPhotoelectricEffect::readQabs() const
{
	///////////////////// Begin copy-paste from SKIRT
	bool reverse = true;
	bool skip1 = false;
	bool skip2 = false;
	bool skip3 = false;
	size_t _Nlambda, _Na;
	// open the file
	ifstream file;
	if (_carbonaceous)
		file = IOTools::ifstreamRepoFile("dat/Gra_81.dat");
	else
		file = IOTools::ifstreamRepoFile("dat/suvSil_81.dat");
	if (!file)
		cout << "Grain data not found!" << endl;

	// skip header lines and read the grid size
	string line;
	while (file.peek() == '#')
		getline(file, line);
	file >> _Na;
	getline(file, line); // ignore anything else on this line
	file >> _Nlambda;
	getline(file, line); // ignore anything else on this line

	// resize the vectors
	_filelambdav.resize(_Nlambda);
	_fileav.resize(_Na);
	_Qabsvv.resize(_Nlambda, vector<double>(_Na));
	_Qscavv.resize(_Nlambda, vector<double>(_Na));
	_asymmparvv.resize(_Nlambda, vector<double>(_Na));

	// determine the loop conditions for wavelength lines
	int kbeg = reverse ? _Nlambda - 1 : 0;
	int kend = reverse ? -1 : _Nlambda;
	int kinc = reverse ? -1 : 1;

	// read the data blocks
	double dummy;
	for (size_t i = 0; i < _Na; i++)
	{
		file >> _fileav[i];
		_fileav[i] *= 1e-6; // convert from micron to m
		getline(file, line); // ignore anything else on this line

		for (int k = kbeg; k != kend; k += kinc)
		{
			if (skip1)
				file >> dummy;
			file >> _filelambdav[k];
			_filelambdav[k] *= 1e-6; // convert from micron to m
			if (skip2)
				file >> dummy;
			file >> _Qabsvv[k][i];
			file >> _Qscavv[k][i];
			if (skip3)
				file >> dummy;
			file >> _asymmparvv[k][i];
			getline(file, line); // ignore anything else on this line
		}
	}

	// close the file
	file.close();
	///////////////////// End copy-paste from SKIRT

	// Convert the wavelengths and grain sizes from microns to centimeters
	for (double& d : _filelambdav)
		d *= 100.; // m to cm
	for (double& d : _fileav)
		d *= 100.; // m to cm
}

std::vector<double>
GrainPhotoelectricEffect::generateQabsv(double a, const std::vector<double>& wavelengthv) const
{
	vector<double> Qabs(wavelengthv.size());
	vector<double> QabsFromFileForA(_filelambdav.size());

	// very simple model
	//        for (size_t i = 0; i < wavelength.size(); i++)
	//        {
	//            Qabs[i] = .75;
	//            if (a < 100 * Constant::ANG_CM) Qabs[i] *= a / 100 / Constant::ANG_CM; // this
	//            works pretty well to simulate the leveling off of the heating efficiency at
	//            1000 Ang
	//        }

	if (a <= _fileav[0]) // extrapolate propto a
	{
		for (size_t i = 0; i < _filelambdav.size(); i++)
			QabsFromFileForA[i] = _Qabsvv[i][0] * a / _fileav[0];
	}
	else // interpolated from data
	{
		size_t a_index = TemplatedUtils::index(a, _fileav);
		double normalDistance = (a - _fileav[a_index - 1]) /
		                        (_fileav[a_index] - _fileav[a_index - 1]);
		// interpolate the values from the file for a specific grain size
		for (size_t i = 0; i < _filelambdav.size(); i++)
			QabsFromFileForA[i] = _Qabsvv[i][a_index - 1] * (1 - normalDistance) +
			                      _Qabsvv[i][a_index] * normalDistance;
	}

#ifdef EXACTGRID
	return QabsFromFileForA;
#endif

	Qabs = TemplatedUtils::linearResample<vector<double>>(QabsFromFileForA, _filelambdav,
	                                                      wavelengthv, -1, -1);
#ifdef PLOT_QABS
	std::stringstream filename;
	filename << "photoelectric/multi-qabs/qabs_a" << setfill('0') << setw(8) << setprecision(2)
	         << fixed << a / Constant::ANG_CM << ".txt";
	ofstream qabsfile = IOTools::ofstreamFile(filename.str());
	for (size_t i = 0; i < wavelengthv.size(); i++)
		qabsfile << wavelengthv[i] * Constant::CM_UM << '\t' << Qabs[i] << endl;
	qabsfile.close();
#endif
	return Qabs;
}

double GrainPhotoelectricEffect::ionizationPotential(double a, int Z) const
{
	double e2a = Constant::ESQUARE / a;
	double ip_v = (Z + .5) * e2a;
	if (Z >= 0)
	{
		// use the same expression for carbonaceous and silicate
		ip_v += _workFunction + (Z + 2) * e2a * 0.3 * Constant::ANG_CM / a; // WD01 eq 2
	}
	else if (_carbonaceous)
	{
		ip_v += _workFunction - e2a * 4.e-8 / (a + 7 * Constant::ANG_CM); // WD01 eq 4
	}
	else // if silicate
	{
		ip_v += 3 / Constant::ERG_EV; // WD01 eq 5
	}
	return ip_v;
}

void GrainPhotoelectricEffect::chargeBalance(double a, const Environment& env,
                                               const std::vector<double>& Qabs, int& resultZmax,
                                               int& resultZmin, std::vector<double>& resultfZ) const
{
	// Express a in angstroms
	double aA = a / Constant::ANG_CM;

	// Shortest wavelength = highest possible energy of a photon
	double hnumax = Constant::PLANCKLIGHT / env.wavelengthv[0];

	/* The maximum charge is one more than the highest charge which still allows ionization by
	   photons of hnumax. */
	resultZmax = floor(
	                ((hnumax - _workFunction) * Constant::ERG_EV / 14.4 * aA + .5 - .3 / aA) /
	                (1 + .3 / aA));

	// The minimum charge is the most negative charge for which autoionization does not occur
	resultZmin = minimumCharge(a);

	if (resultZmax < resultZmin)
		Error::runtime("Zmax is smaller than Zmin");
	resultfZ.resize(resultZmax - resultZmin + 1, -1);

	/* We will cut off the distribution at some point past the maximum (in either the positive
	   or the negative direction). */
	double lowerLimit = 1.e-6;
	int trimLow = 0;
	int trimHigh = resultfZ.size() - 1;

	int centerZ = 0;

	// Find the central peak using a binary search
	int upperbound = resultZmax;
	int lowerbound = resultZmin;
	int current = floor((resultZmax + resultZmin) / 2.);
	while (current != lowerbound)
	{
		// Upward ratio in detailed balance equation
		double Jpe = emissionRate(a, current, env.wavelengthv, Qabs, env.energyDensityv);
		double Jion = collisionalChargingRate(a, env.T, current, 1, Constant::PROTONMASS,
		                                      env.np);
		double Je = collisionalChargingRate(a, env.T, current + 1, -1,
		                                    Constant::ELECTRONMASS, env.ne);
		double factor = (Jpe + Jion) / Je;

		// Upward slope => the maximum is more to the right
		// Downward slope => the maximum is more to the left
		if (factor > 1)
			lowerbound = current;
		else
			upperbound = current;
		// Move the cursor to the center of the new bounds
		current = floor((lowerbound + upperbound) / 2.);
		// cout << "Searching maximum in ("<<lowerbound<<" | "<<current<<" |
		// "<<upperbound<<")"<<endl;
	}
	centerZ = current;

	// Apply detailed balance equation ...
	bool passedMaximum = false;
	double maximum = 0;
	resultfZ[centerZ - resultZmin] = 1;
	// ... for Z > centerZ
	for (int z = centerZ + 1; z <= resultZmax; z++)
	{
		double Jpe = emissionRate(a, z - 1, env.wavelengthv, Qabs, env.energyDensityv);
		double Jion = collisionalChargingRate(a, env.T, z - 1, 1, Constant::PROTONMASS,
		                                      env.np);
		double Je = collisionalChargingRate(a, env.T, z, -1, Constant::ELECTRONMASS,
		                                    env.ne);
		// cout << "z = "<<z<<" Jpe = "<<Jpe<<" Jion = "<<Jion<<" Je = "<<Je<<endl;
		int index = z - resultZmin;
		resultfZ[index] = resultfZ[index - 1] * (Jpe + Jion) / Je;

		// Cut off calculation once the values past the maximum have a relative size of 1e-6
		// or less Detects the maximum
		if (!passedMaximum && resultfZ[index] < resultfZ[index - 1])
		{
			passedMaximum = true;
			maximum = max(maximum, resultfZ[index - 1]);
		}
		// Detects the cutoff point once the maximum has been found
		if (passedMaximum && resultfZ[index] / maximum < lowerLimit)
		{
			trimHigh = index;
			break;
		}
	}

	// ... for Z < centerZ
	for (int z = centerZ - 1; z >= resultZmin; z--)
	{
		double Jpe = emissionRate(a, z, env.wavelengthv, Qabs, env.energyDensityv);
		double Jion = collisionalChargingRate(a, env.T, z, 1, Constant::PROTONMASS, env.np);
		double Je = collisionalChargingRate(a, env.T, z + 1, -1, Constant::ELECTRONMASS,
		                                    env.np);
		int index = z - resultZmin;
		resultfZ[index] = resultfZ[index + 1] * Je / (Jpe + Jion);

		if (!passedMaximum && resultfZ[index] < resultfZ[index + 1])
		{
			passedMaximum = true;
			maximum = max(maximum, resultfZ[index + 1]);
		}
		if (passedMaximum && resultfZ[index] / maximum < lowerLimit)
		{
			trimLow = index;
			break;
		}
	}

	// Trim the top first! Makes things (i.e. dealing with the indices) much easier.
	resultfZ.erase(resultfZ.begin() + trimHigh + 1, resultfZ.end());
	resultZmax = resultZmin + trimHigh;

	resultfZ.erase(resultfZ.begin(), resultfZ.begin() + trimLow);
	resultZmin += trimLow;

	// Normalize
	double sum = 0.;
	for (double d : resultfZ)
		sum += d;
	for (double& d : resultfZ)
		d /= sum;

	// Test the result
	//    for (int z = resultZmin + 1; z <= resultZmax - 1; z++){
	//        double Jpe = emissionRate(a, z, wavelength, Qabs, isrf);
	//        double Jion = collisionalChargingRate(a, z, 1, Constant::HMASS_CGS);
	//        double Je = collisionalChargingRate(a, z, -1, Constant::ELECTRONMASS);
	//        int index = z - resultZmin;
	//        double departure = resultfZ[index] * (Jpe + Jion + Je);

	//        Jpe = emissionRate(a, z-1, wavelength, Qabs, isrf);
	//        Jion = collisionalChargingRate(a, z-1, 1, Constant::HMASS_CGS);
	//        double arrival = resultfZ[index-1] * (Jpe + Jion);

	//        Je = collisionalChargingRate(a, z+1, -1, Constant::ELECTRONMASS);
	//        arrival += resultfZ[index+1] * Je;

	//        cout << "Total change rate in population z = "<<z<<" : "<<arrival -
	//        departure<<endl;
	//    }
}

double GrainPhotoelectricEffect::heatingRateAZ(double a, int Z, const Array& wavelengthv,
                                                 const std::vector<double>& Qabs,
                                                 const Array& energyDensity_lambda) const
{
	const double e2a = Constant::ESQUARE / a;
	size_t nLambda = wavelengthv.size();

	// Two separate integrands for photoelectric effect and photodetachment
	vector<double> peIntegrandv(nLambda, 0);
	vector<double> pdIntegrandv(nLambda, 0);

	// Quantities independent of nu
	double ip_v = ionizationPotential(a, Z);

	double Emin = Z >= 0 ? 0
	                     : -(Z + 1) * e2a / (1 + std::pow(27. * Constant::ANG_CM / a,
	                                                      0.75)); // WD01 eq 7

	double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin; // WD01 eq 6

	double Elow = Z < 0 ? Emin : -(Z + 1) * e2a; // WD01 text between eq 10 and 11

	for (size_t lambda_index = 0; lambda_index < nLambda; lambda_index++)
	{
		double hnu = Constant::PLANCKLIGHT / wavelengthv[lambda_index];

		// No contribution below the photoelectric threshold
		if (hnu > hnu_pet)
		{
			// Quantities dependent on nu
			double hnuDiff = hnu - hnu_pet;

			double Emax = hnuDiff + Emin;
			double Ehigh = Z < 0 ? Emax : hnuDiff;

			// Calculate y2 from eq 11
			double Ediff = Ehigh - Elow;
			double y2 = Z >= 0
			                            ? Ehigh * Ehigh * (Ehigh - 3 * Elow) / Ediff /
			                                              Ediff / Ediff
			                            : 1;

			// The integral over the electron energy distribution
			double IntE = energyIntegral(Elow, Ehigh, Emin, Emax);
			double Y = yield(a, Z, hnuDiff, Elow, Ehigh);

			peIntegrandv[lambda_index] = Y * Qabs[lambda_index] *
			                             energyDensity_lambda[lambda_index] / hnu *
			                             IntE / y2;
		}

		// If applicable, also calculate integrand for photodetachment
		if (Z < 0)
		{
			// Photodetachment
			double hnu_pdt = ip_v + Emin; // eq 18

			// No contribution below the photodetachment threshold
			if (hnu > hnu_pdt)
			{
				double DeltaE = 3 / Constant::ERG_EV;
				double x = (hnu - hnu_pdt) / DeltaE;
				double denom = 1 + x * x / 3;
				double sigma_pdt = x / denom / denom; // eq 20

				pdIntegrandv[lambda_index] = sigma_pdt *
				                             energyDensity_lambda[lambda_index] /
				                             hnu * (hnu - hnu_pdt + Emin);
			}
		}
	}

	double peIntegral = Constant::PI * a * a * Constant::LIGHT *
	                    TemplatedUtils::integrate<double>(wavelengthv, peIntegrandv);

	// Constant factor from sigma_pdt moved in front of integral
	double pdIntegral =
	                Z < 0
	                                ? 1.2e-17 * (-Z) * Constant::LIGHT *
	                                                  TemplatedUtils::integrate<double>(
	                                                                  wavelengthv, pdIntegrandv)
	                                : 0;

	return peIntegral + pdIntegral;
}

double GrainPhotoelectricEffect::heatingRateA(double a, const Environment& env,
                                                const std::vector<double>& Qabs) const
{
	double totalHeatingForGrainSize = 0;

	vector<double> fZ;
	int Zmin, Zmax;
	chargeBalance(a, env, Qabs, Zmax, Zmin, fZ);

	printf("Z in (%d, %d)\n", Zmin, Zmax);

	for (int Z = Zmin; Z <= Zmax; Z++)
	{
		double fZz = fZ[Z - Zmin];
		double heatAZ = heatingRateAZ(a, Z, env.wavelengthv, Qabs, env.energyDensityv);
		totalHeatingForGrainSize += fZz * heatAZ;
	}

	// The net heating rate (eq 41 without denominator)
	totalHeatingForGrainSize -= recombinationCoolingRate(a, env, fZ, Zmin);

#ifdef PLOT_FZ
	std::stringstream filename;
	filename << "photoelectric/multi-fz/fz_a" << setfill('0') << setw(8) << setprecision(2)
	         << fixed << a / Constant::ANG_CM << ".txt";
	ofstream outvar = IOTools::ofstreamFile(filename.str());
	outvar << "# carbon = " << _carbonaceous << endl;
	outvar << "# a = " << a << endl;
	outvar << "# Teff = " << _Tc << endl;
	outvar << "# G0 = " << _G0 << endl;
	outvar << "# ne = " << _electronDensity << endl;
	outvar << "# Tgas = " << _gasTemperature << endl;

	for (int z = Zmin; z <= Zmax; z++)
		outvar << z << '\t' << fZ[z - Zmin] << '\n';
	outvar.close();
#endif

	return totalHeatingForGrainSize;
}

double GrainPhotoelectricEffect::heatingRate(const Environment& env,
                                               const std::vector<double>& grainSizev,
                                               const std::vector<double>& grainDensityv,
                                               const std::vector<std::vector<double>>& absQvv) const
{
	double total;
	for (size_t m = 0; m < grainSizev.size(); m++)
	{
		total += grainDensityv[m] * heatingRateA(grainSizev[m], env, absQvv[m]);
	}
	return total;
}

double GrainPhotoelectricEffect::emissionRate(double a, int Z, const Array& wavelengthv,
                                                const std::vector<double>& Qabs,
                                                const Array& energyDensity_lambda) const
{
	const double e2a = Constant::ESQUARE / a;
	size_t nLambda = wavelengthv.size();

	// Calculate the integrandum at the frequencies of the wavelength grid
	vector<double> peIntegrandv(nLambda, 0);
	vector<double> pdIntegrandv(nLambda, 0);

	// Quantities independent of nu
	double ip_v = ionizationPotential(a, Z);

	// WD01 eq 7
	double Emin = Z >= 0 ? 0
	                     : -(Z + 1) * e2a / (1 + std::pow(27. * Constant::ANG_CM / a,
	                                                      0.75));

	// WD01 eq 6
	double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin;

	// WD01 text between eq 10 and 11
	double Elow = Z < 0 ? Emin : -(Z + 1) * e2a;

	for (size_t lambda_index = 0; lambda_index < nLambda; lambda_index++)
	{
		double hnu = Constant::PLANCKLIGHT / wavelengthv[lambda_index];

		// No contribution below the photoelectric threshold
		if (hnu > hnu_pet)
		{
			// Quantities dependent on nu
			double hnuDiff = hnu - hnu_pet;
			double Ehigh = Z < 0 ? Emin + hnuDiff : hnuDiff;

			// The integral over the electron energy distribution
			double Y = yield(a, Z, hnuDiff, Elow, Ehigh);

			peIntegrandv[lambda_index] =
			                Y * Qabs[lambda_index] * energyDensity_lambda[lambda_index] / hnu;
		}

		if (Z < 0)
		{
			// Photodetachment, WD01 eq 18
			double hnu_pdt = ip_v + Emin;

			// No contribution below the photodetachment threshold
			if (hnu > hnu_pdt)
			{
				double DeltaE = 3 / Constant::ERG_EV;
				double x = (hnu - hnu_pdt) / DeltaE;
				double denom = 1 + x * x / 3;
				// Cross section, WD01 eq 20
				double sigma_pdt = x / denom / denom;

				pdIntegrandv[lambda_index] = sigma_pdt * energyDensity_lambda[lambda_index] / hnu;
			}
		}
	}

	double peIntegral = Constant::PI * a * a * Constant::LIGHT *
	                    TemplatedUtils::integrate<double>(wavelengthv, peIntegrandv);

	// Constant factor from sigma_pdt (WD01 eq 20) moved in front of integral
	double pdIntegral =
	                Z < 0
	                                ? 1.2e-17 * (-Z) * Constant::LIGHT *
	                                                  TemplatedUtils::integrate<double>(
	                                                                  wavelengthv, pdIntegrandv)
	                                : 0;
	return peIntegral + pdIntegral;
}

double GrainPhotoelectricEffect::energyIntegral(double Elow, double Ehigh, double Emin,
                                                  double Emax) const
{
	double Ediff = Ehigh - Elow;
	double Ediff3 = Ediff * Ediff * Ediff;

	/* Compute integral f(E)E dE analytically, with f(E) defined by WD01 eq 10 f(E) is a
	   parabola, and therefore f(E)E is a third order polynomial.  Thus the integral of f(E)E dE
	   is a fourth order polynomial: a/4 (max4 - min4) + b/3 (max3 - min3) + c/2 (max2
	   -min2). */
	double Emax2 = Emax * Emax;
	double Emin2 = Emin * Emin;
	return 6 / Ediff3 * (-(Emax2 * Emax2 - Emin2 * Emin2) / 4. +
	                     (Ehigh + Elow) * (Emax2 * Emax - Emin2 * Emin) / 3. -
	                     Elow * Ehigh * (Emax2 - Emin2) / 2.);
}

double GrainPhotoelectricEffect::yield(double a, int Z, double hnuDiff, double Elow,
                                         double Ehigh) const
{
	if (hnuDiff < 0)
		Error::runtime("Frequency is smaller than photoelectric threshold.");

	// Compute yield (y2, y1, y0, and finally Y)

	// Calculate y2 from eq 11
	double Ediff = Ehigh - Elow;
	double y2 = Z >= 0 ? Ehigh * Ehigh * (Ehigh - 3. * Elow) / Ediff / Ediff / Ediff : 1;

	// Calculate y1 from grain properties and eq 13, 14
	// double imaginaryRefIndex = 1; // should be wavelength-dependent
	// double la = Constant::LIGHT / nu / 4 / Constant::PI / imaginaryRefIndex;
	const double la = 100 * Constant::ANG_CM; // value from 1994-Bakes. WD01 uses the above one
	double beta = a / la;
	const double le = 10 * Constant::ANG_CM;
	double alpha = beta + a / le;

	double alpha2 = alpha * alpha;
	double beta2 = beta * beta;
	double y1 = beta2 / alpha2 * (alpha2 - 2. * alpha - 2. * expm1(-alpha)) /
	            (beta2 - 2. * beta - 2. * expm1(-beta));

	// Calculate y0 from eq 9, 16, 17
	double thetaOverW = hnuDiff;
	if (Z >= 0)
		thetaOverW += (Z + 1) * Constant::ESQUARE / a;
	thetaOverW /= _workFunction;
	double thetaOverWtothe5th = thetaOverW * thetaOverW;
	thetaOverWtothe5th *= thetaOverWtothe5th * thetaOverW;

	double y0;
	// Different formulae for carbonaceous and silicate
	if (_carbonaceous)
		y0 = 9e-3 * thetaOverWtothe5th / (1 + 3.7e-2 * thetaOverWtothe5th);
	else
		y0 = 0.5 * thetaOverW / (1. + 5. * thetaOverW);

	return y2 * min(y0 * y1, 1.);
}

int GrainPhotoelectricEffect::minimumCharge(double a) const
{
	double aA = a / Constant::ANG_CM;

	// WD01 eq 23
	double Uait = _carbonaceous ? 3.9 + 0.12 * aA + 2. / aA : 2.5 + 0.07 * aA + 8. / aA;
	// WD01 eq 24
	return floor(-Uait / 14.4 * aA) + 1;
}

double GrainPhotoelectricEffect::stickingCoefficient(double a, int Z, int z_i) const
{
	// ions
	if (z_i >= 0)
		return 1;
	// electrons
	else if (z_i == -1)
	{
		// electron mean free path length in grain
		double le = 10. * Constant::ANG_CM;
		double pElasticScatter = .5;
		// negative and neutral grains
		if (Z <= 0)
		{
			if (Z > minimumCharge(a))
			{
				// number of carbon atoms
				double NC = 468 * a * a * a / 1.e-21;
				return 0.5 * (-expm1(-a / le)) / (1 + exp(20 - NC));
			}
			else
				return 0;
		}
		// positive grains
		else
			return (1 - pElasticScatter) * (-expm1(-a / le));
	}
	else
		return 0;
}

double GrainPhotoelectricEffect::collisionalChargingRate(double a, double gasT, int Z,
                                                           int particleCharge, double particleMass,
                                                           double particleDensity) const
{
	double stick = stickingCoefficient(a, Z, particleCharge);

	double kT = Constant::BOLTZMAN * gasT;

	// WD eq 26: akT / q^2 = akT / e^2 / z^2
	double tau = a * kT / Constant::ESQUARE / particleCharge /
	             particleCharge;
	// Ze / q = Z / z
	double ksi = static_cast<double>(Z) /
	             static_cast<double>(particleCharge);

	double Jtilde;
	if (ksi < 0)
	{
		Jtilde = (1. - ksi / tau) * (1. + sqrt(2. / (tau - 2. * ksi)));
	}
	else if (ksi > 0)
	{
		double toSquare = 1. + 1. / sqrt(4. * tau + 3. * ksi);
		Jtilde = toSquare * toSquare * exp(-ksi / (1. + 1. / sqrt(ksi)) / tau);
	}
	else
	{
		Jtilde = 1. + sqrt(Constant::PI / 2. / tau);
	}
	return particleDensity * stick * sqrt(8. * kT * Constant::PI / particleMass) * a * a *
	       Jtilde;
}

double GrainPhotoelectricEffect::lambdaTilde(double tau, double ksi) const
{
	// Found in 1987-Draine-Sutin
	if (ksi < 0)
	{
		return (2. - ksi / tau) * (1. + 1. / sqrt(tau - ksi));
	}
	else if (ksi > 0)
	{
		return (2. + ksi / tau) * (1. + 1. / sqrt(3. / 2. / tau + 3. * ksi)) *
		       exp(-ksi / (1. + 1. / sqrt(ksi)) / tau);
	}
	else
	{
		return 2. + 3. / 2. * sqrt(Constant::PI / 2. / tau);
	}
}

double GrainPhotoelectricEffect::recombinationCoolingRate(double a, const Environment& env,
                                                            const std::vector<double>& fZ,
                                                            int Zmin) const
{
	// Calculates WD01 equation 42
	double kT = Constant::BOLTZMAN * env.T;
	double eightkT3DivPi = 8 * kT * kT * kT / Constant::PI;

	int Zmax = Zmin + fZ.size() - 1;

	// For every collision partner, add the contibutions for each possible grain charge.
	double particleSum = 0;
	for (size_t i = 0; i < env.particleChargev.size(); i++)
	{
		// tau = akT / q^2 (WD01 eq 26)
		int z_i = env.particleChargev[i];
		double tau = a * kT / z_i / z_i / Constant::ESQUARE;

		double Zsum = 0;
		for (int Z = Zmin; Z <= Zmax; Z++)
		{
			// ksi = Ze / q_i = Z / z_i
			double ksi = Z / static_cast<double>(z_i);
			Zsum += stickingCoefficient(a, Z, z_i) * fZ[Z - Zmin] *
			        lambdaTilde(tau, ksi);
		}
		particleSum += env.particleDensityv[i] *
		               sqrt(eightkT3DivPi / env.particleMassv[i]) * Zsum;
	}

	// Previous implementation, which was correct. Kept for reference at the moment.
	//	// electrons
	//	double Zsum = 0;
	//	for (int z = Zmin; z <= Zmax; z++)
	//	{
	//		double ksi = -z; // ksi = Ze / q_e = -Z
	//		Zsum += stickingCoefficient(a, z, -1) * fZ[z - Zmin] * lambdaTilde(tau, ksi);
	//	}
	//	particleSum += _electronDensity * sqrt(eightkT3DivPi / Constant::ELECTRONMASS) * Zsum;
	//
	//	// (H+) ions
	//	Zsum = 0;
	//	for (int z = Zmin; z <= Zmax; z++)
	//	{
	//		double ksi = z; // ksi = Ze / q_i = z
	//		Zsum += stickingCoefficient(a, z, 1) * fZ[z - Zmin] * lambdaTilde(tau, ksi);
	//	}
	//	particleSum += _electronDensity * sqrt(eightkT3DivPi / Constant::PROTONMASS) * Zsum;

	/* The second term of equation 42: autoionization of grains with the most negative charge
	   inhibits the cooling of the gas. */
	// EA(Zmin) = IP(Zmin-1) because IP(Z) = EA(Z+1)
	double secondTerm = 0;
	/* This term is only included when the population of the maximally negative grain charge
	   minimumCharge is significant. If it is not siginicant, then fZ will not cover
	   minimumCharge, (and Zmin > minimumCharge). */
	if (Zmin == minimumCharge(a))
		double secondTerm = fZ[0] *
		                    collisionalChargingRate(a, env.T, Zmin, -1,
		                                            Constant::ELECTRONMASS, env.ne) *
		                    ionizationPotential(a, Zmin - 1);

	return Constant::PI * a * a * particleSum + secondTerm;
}

double GrainPhotoelectricEffect::yieldFunctionTest() const
{
	// Parameters
	const int Z = 10;
	vector<double> av = {4e-8, 10e-8, 30e-8, 100e-8, 300e-8};

	// Plot range
	const double hnuMin = 5 / Constant::ERG_EV;
	const double hnuMax = 15 / Constant::ERG_EV;
	const size_t N = 500;

	ofstream out = IOTools::ofstreamFile("photoelectric/yieldTest.dat");
	for (double a : av)
	{
		out << "# a = " << a << '\n';
		const double e2a = Constant::ESQUARE / a;

		// Quantities independent of nu
		double ip_v = ionizationPotential(a, Z);

		double Emin = Z >= 0 ? 0
		                     : -(Z + 1) * e2a / (1 + std::pow(27. * Constant::ANG_CM / a,
		                                                      0.75)); // WD01 eq 7

		double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin; // WD01 eq 6

		double Elow = Z < 0 ? Emin : -(Z + 1) * e2a; // WD01 text between eq 10 and 11

		double hnu = hnuMin;
		const double step = (hnuMax - hnuMin) / N;
		for (size_t n = 0; n < N; n++)
		{
			double hnuDiff = hnu - hnu_pet;
			double Ehigh = Z < 0 ? hnuDiff + Emin : hnuDiff;

			if (hnuDiff > 0)
				out << hnu * Constant::ERG_EV << '\t'
				    << yield(a, Z, hnuDiff, Elow, Ehigh) << '\n';
			hnu += step;
		}
		out << '\n';
	}
	out.close();
	return 0.0;
}

double GrainPhotoelectricEffect::heatingRateTest(double G0, double gasT, double ne) const
{
	// Wavelength grid
	const vector<double>& frequencyv = Testing::generateGeometricGridv(
	                _nWav, Constant::LIGHT / _maxWav, Constant::LIGHT / _minWav);

	// Input spectrum
	const Array& specificIntensityv = Testing::generateSpecificIntensityv(frequencyv, _Tc, G0);

	// Convert to wavelength units
	const vector<double>& wavelengthv = Testing::freqToWavGrid(frequencyv);
	Array energyDensity_lambda =
	                Constant::FPI / Constant::LIGHT *
	                Testing::freqToWavSpecificIntensity(frequencyv, specificIntensityv);

	// Write out wavelengths and isrf
	ofstream isrfOf = IOTools::ofstreamFile("photoelectric/isrf_ulambda.dat");
	for (size_t i = 0; i < wavelengthv.size(); i++)
		isrfOf << wavelengthv[i] * Constant::CM_UM << '\t' << energyDensity_lambda[i]
		       << endl;
	isrfOf.close();

	// Gather environment parameters
	const Environment env(Array(wavelengthv.data(), wavelengthv.size()), energyDensity_lambda,
	                      gasT, ne, ne, {-1, 1}, {ne, ne},
	                      {Constant::ELECTRONMASS, Constant::PROTONMASS});

	// Read absorption efficiency from SKIRT file into local memory
	readQabs();

	/* File that writes out the absorption efficiency, averaged using the input radiation field
	   as weights. */
	ofstream avgQabsOf = IOTools::ofstreamFile("photoelectric/avgQabsInterp.txt");

	// Grain sizes for test
	double aMin = 1.5 * Constant::ANG_CM;
	double aMax = 10000 * Constant::ANG_CM;
	const size_t Na = 60;
	double aStepFactor = std::pow(aMax / aMin, 1. / Na);
	double a = aMin;

	// Output file will contain one line for every grain size
	stringstream efficiencyFnSs;
	efficiencyFnSs << "photoelectric/efficiencyG" << setprecision(4) << scientific << G0
	               << ".dat";
	std::ofstream efficiencyOf = IOTools::ofstreamFile(efficiencyFnSs.str());

	// For every grain size
	for (size_t m = 0; m < Na; m++)
	{
		// From the data read in from the SKIRT file, interpolate an absorption efficiency.
		vector<double> Qabs = generateQabsv(a, wavelengthv);

		// Integrate over the radiation field
		vector<double> ulambdaTimesQabs(Qabs.size());
		for (size_t n = 0; n < Qabs.size(); n++)
			ulambdaTimesQabs[n] = Qabs[n] * energyDensity_lambda[n];
		double uTimesQabsIntegral =
		                TemplatedUtils::integrate<double>(wavelengthv, ulambdaTimesQabs);

		cout << "Size " << a / Constant::ANG_CM << endl;

		// Calculate and write out the heating efficiency
		double heating = GrainPhotoelectricEffect::heatingRateA(a, env, Qabs);
		double totalAbsorbed = Constant::PI * a * a * Constant::LIGHT * uTimesQabsIntegral;
		efficiencyOf << a / Constant::ANG_CM << '\t' << heating / totalAbsorbed << '\n';

		// Calculate and write out the ISRF-averaged absorption efficiency
		double uIntegral = TemplatedUtils::integrate<double>(wavelengthv,
		                                                     energyDensity_lambda);
		double avgQabs = uTimesQabsIntegral / uIntegral;
		avgQabsOf << a / Constant::ANG_CM << '\t' << avgQabs << endl;
		a *= aStepFactor;
	}
	efficiencyOf.close();
	cout << "Wrote " << efficiencyFnSs.str() << endl;
	avgQabsOf.close();
	cout << "Wrote avgQabsInterp.txt" << endl;
	cout << "Charging parameter = " << G0 * sqrt(gasT) / ne << endl;
	return 0.0;
}

double GrainPhotoelectricEffect::chargeBalanceTest(double G0, double gasT, double ne,
                                                     double np) const
{
	readQabs();

// Wavelength grid
#ifdef EXACTGRID
	vector<double> wavelengthv = _filelambdav;
#else
	vector<double> wavelengthv = Testing::generateGeometricGridv(_nWav, _minWav, _maxWav);
#endif

	// Input spectrum
	Array isrfv = Testing::generateSpecificIntensityv(wavelengthv, _Tc, G0);

	const Environment env(Array(wavelengthv.data(), wavelengthv.size()), isrfv, gasT, ne, np,
	                      {-1, 1}, {ne, np}, {Constant::ELECTRONMASS, Constant::PROTONMASS});

	// Grain size
	double a = 200. * Constant::ANG_CM;

	// Absorption spectrum
	vector<double> Qabsv = generateQabsv(a, wavelengthv);

	// Calculate charge distribution
	vector<double> fZv;
	int Zmax, Zmin;
	chargeBalance(a, env, Qabsv, Zmax, Zmin, fZv);

	std::cout << "Zmax = " << Zmax << " Zmin = " << Zmin << " len fZ = " << fZv.size()
	          << std::endl;

	std::ofstream out = IOTools::ofstreamFile("photoelectric/fZ.txt");
	out << "# carbon = " << _carbonaceous << endl;
	out << "# a = " << a << endl;
	out << "# Teff = " << _Tc << endl;
	out << "# G0 = " << G0 << endl;
	out << "# ne = " << ne << endl;
	out << "# Tgas = " << gasT << endl;

	for (int z = Zmin; z <= Zmax; z++)
		out << z << '\t' << fZv[z - Zmin] << '\n';
	out.close();

	return 0.;
}
