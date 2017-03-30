#include "PhotoelectricHeating.h"

#include "NumUtils.h"
#include "Testing.h"
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

//#define EXACTGRID

using namespace std;

namespace
{
vector<double> _filelambdav, _fileav;
vector<vector<double>> _Qabsvv, _Qscavv, _asymmparvv;
}
void PhotoelectricHeatingRecipe::readQabs() const
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
		file.open("../git/dat/Gra_81.dat");
	else
		file.open("../git/dat/suvSil_81.dat");
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
		_fileav[i] *= 1e-6;  // convert from micron to m
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
PhotoelectricHeatingRecipe::generateQabsv(double a, const std::vector<double>& wavelengthv) const
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
		size_t a_index = NumUtils::index<double>(a, _fileav);
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

	Qabs = NumUtils::interpol<double>(QabsFromFileForA, _filelambdav, wavelengthv, -1, -1);
#ifdef PLOT_QABS
	ofstream qabsfile;
	std::stringstream filename;
	filename << "/Users/drvdputt/GasModule/run/multi-qabs/qabs_a" << setfill('0') << setw(8)
	         << setprecision(2) << fixed << a / Constant::ANG_CM << ".txt";
	qabsfile.open(filename.str());
	for (size_t i = 0; i < wavelengthv.size(); i++)
		qabsfile << wavelengthv[i] * Constant::CM_UM << '\t' << Qabs[i] << endl;
	qabsfile.close();
#endif
	return Qabs;
}

double PhotoelectricHeatingRecipe::ionizationPotential(double a, int Z) const
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

double
PhotoelectricHeatingRecipe::heatingRateAZ(double a, int Z, const std::vector<double>& wavelengthv,
                                          const std::vector<double>& Qabs,
                                          const std::vector<double>& energyDensity_lambda) const
{
	const double e2a = Constant::ESQUARE / a;
	size_t nLambda = wavelengthv.size();

	// Two separate integrands for photoelectric effect and photodetachment
	vector<double> peIntegrandv(nLambda, 0);
	vector<double> pdIntegrandv(nLambda, 0);

	// Quantities independent of nu
	double ip_v = ionizationPotential(a, Z);

	double Emin = Z >= 0 ? 0
	                     : -(Z + 1) * e2a /
	                                              (1 + std::pow(27. * Constant::ANG_CM / a,
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
			double y2 = Z >= 0 ? Ehigh * Ehigh * (Ehigh - 3 * Elow) / Ediff / Ediff /
			                                            Ediff
			                   : 1;

			// The integral over the electron energy distribution
			double IntE = energyIntegral(Elow, Ehigh, Emin, Emax);
			double Y = yield(a, Z, hnuDiff, Elow, Ehigh);

			peIntegrandv[lambda_index] = Y * Qabs[lambda_index] *
			                             energyDensity_lambda[lambda_index] / hnu *
			                             IntE / y2;
		}

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
	                    NumUtils::integrate<double>(wavelengthv, peIntegrandv);

	// Constant factor from sigma_pdt moved in front of integral
	double pdIntegral =
	                Z < 0 ? 1.2e-17 * (-Z) * Constant::LIGHT *
	                                                NumUtils::integrate<double>(wavelengthv,
	                                                                            pdIntegrandv)
	                      : 0;

	return peIntegral + pdIntegral;
}

double
PhotoelectricHeatingRecipe::heatingRateA(double a, const std::vector<double>& wavelengthv,
                                         const std::vector<double>& Qabs,
                                         const std::vector<double>& energyDensity_lambda) const
{
	double totalHeatingForGrainSize = 0;

	vector<double> fZ;
	int Zmin, Zmax;
	chargeBalance(a, wavelengthv, Qabs, energyDensity_lambda, Zmax, Zmin, fZ);

	printf("Z in (%d, %d)\n", Zmin, Zmax);

	for (int z = Zmin; z <= Zmax; z++)
	{
		double fZz = fZ[z - Zmin];
		double heatAZ = heatingRateAZ(a, z, wavelengthv, Qabs, energyDensity_lambda);
		totalHeatingForGrainSize += fZz * heatAZ;
	}

	// The net heating rate
	totalHeatingForGrainSize -=
	                recombinationCoolingRate(a, fZ, Zmin); // eq 41 without denominator

#ifdef PLOT_FZ
	std::ofstream outvar;
	std::stringstream filename;
	filename << "/Users/drvdputt/GasModule/run/multi-fz/fz_a" << setfill('0') << setw(8)
	         << setprecision(2) << fixed << a / Constant::ANG_CM << ".txt";
	outvar.open(filename.str());
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

double PhotoelectricHeatingRecipe::heatingRate(const std::vector<double>& wavelengthv,
                                               const std::vector<double>& Qabs,
                                               const std::vector<double>& isrf) const
{
	// Need size distribution

	return 0.;
}

double PhotoelectricHeatingRecipe::emissionRate(double a, int Z,
                                                const std::vector<double>& wavelengthv,
                                                const std::vector<double>& Qabs,
                                                const std::vector<double>& isrf) const
{
	const double e2a = Constant::ESQUARE / a;
	size_t nLambda = wavelengthv.size();

	// Calculate the integrandum at the frequencies of the wavelength grid
	vector<double> peIntegrandv(nLambda, 0);
	vector<double> pdIntegrandv(nLambda, 0);

	// Quantities independent of nu
	double ip_v = ionizationPotential(a, Z);

	double Emin = Z >= 0 ? 0
	                     : -(Z + 1) * e2a /
	                                              (1 + std::pow(27. * Constant::ANG_CM / a,
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
			double Ehigh = Z < 0 ? Emin + hnuDiff : hnuDiff;

			// The integral over the electron energy distribution
			double Y = yield(a, Z, hnuDiff, Elow, Ehigh);

			peIntegrandv[lambda_index] =
			                Y * Qabs[lambda_index] * isrf[lambda_index] / hnu;
		}

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
				double sigma_pdt = x / denom / denom; // eq 20x

				pdIntegrandv[lambda_index] = sigma_pdt * isrf[lambda_index] / hnu;
			}
		}
	}

	double peIntegral = Constant::PI * a * a * Constant::LIGHT *
	                    NumUtils::integrate<double>(wavelengthv, peIntegrandv);

	// Constant factor from sigma_pdt moved in front of integral
	double pdIntegral =
	                Z < 0 ? 1.2e-17 * (-Z) * Constant::LIGHT *
	                                                NumUtils::integrate<double>(wavelengthv,
	                                                                            pdIntegrandv)
	                      : 0;
	return peIntegral + pdIntegral;
}

double PhotoelectricHeatingRecipe::energyIntegral(double Elow, double Ehigh, double Emin,
                                                  double Emax) const
{
	double Ediff = Ehigh - Elow;
	double Ediff3 = Ediff * Ediff * Ediff;

	// Compute integral f(E)E dE analytically, where f(E) is the parabola defined by WD01 eq 10
	// integral f(E)dE is of the form a/4 (max4 - min4) + b/3 (max3 - min3) + c/2 (max2 -min2)
	double Emax2 = Emax * Emax;
	double Emin2 = Emin * Emin;
	return 6 / Ediff3 *
	       (-(Emax2 * Emax2 - Emin2 * Emin2) / 4. +
	        (Ehigh + Elow) * (Emax2 * Emax - Emin2 * Emin) / 3. -
	        Elow * Ehigh * (Emax2 - Emin2) / 2.);
}

double PhotoelectricHeatingRecipe::yield(double a, int Z, double hnuDiff, double Elow,
                                         double Ehigh) const
{
	if (hnuDiff < 0)
		throw std::range_error("Frequency is smaller than photoelectric threshold.");

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

int PhotoelectricHeatingRecipe::minimumCharge(double a) const
{
	double aA = a / Constant::ANG_CM;

	// WD01 eq 23
	double Uait = _carbonaceous ? 3.9 + 0.12 * aA + 2. / aA : 2.5 + 0.07 * aA + 8. / aA;
	// WD01 eq 24
	return floor(-Uait / 14.4 * aA) + 1;
}

double PhotoelectricHeatingRecipe::stickingCoefficient(double a, int Z, int z_i) const
{
	// electrons
	if (z_i == -1)
	{
		double le = 10. * Constant::ANG_CM;
		double pElasticScatter = .5;
		// negative and neutral grains
		if (Z <= 0)
		{
			if (Z > minimumCharge(a))
			{
				// electron mean free path length in grain
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
	// ions
	else
		return 1;
}

double PhotoelectricHeatingRecipe::collisionalChargingRate(double a, int Z, int z_i,
                                                           double m_i) const
{
	double stick = stickingCoefficient(a, Z, z_i);

	double kT = Constant::BOLTZMAN * _gasTemperature;

	double tau = a * kT / Constant::ESQUARE / z_i /
	             z_i; // notes eq 38 (akT / q = akT / e^2 / z^2)
	double ksi = static_cast<double>(Z) /
	             static_cast<double>(z_i); // notes eq 39 (Ze / q = Z / z)

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
	return _electronDensity * stick * sqrt(8. * kT * Constant::PI / m_i) * a * a * Jtilde;
}

double PhotoelectricHeatingRecipe::lambdaTilde(double tau, double ksi) const
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

double PhotoelectricHeatingRecipe::recombinationCoolingRate(double a, const std::vector<double>& fZ,
                                                            int Zmin) const
{
	double kT = Constant::BOLTZMAN * _gasTemperature;
	double eightkT3DivPi = 8 * kT * kT * kT / Constant::PI;

	int Zmax = Zmin + fZ.size() - 1;

	double particleSum = 0;

	// tau = akT / q^2 = akT / e^2 (this needs to change when ions with charges > 1 are
	// introduced)
	double tau = a * kT / Constant::ESQUARE;

	// electrons
	double Zsum = 0;
	for (int z = Zmin; z <= Zmax; z++)
	{
		double ksi = -z; // ksi = Ze / q_e = -Z
		Zsum += stickingCoefficient(a, z, -1) * fZ[z - Zmin] * lambdaTilde(tau, ksi);
	}
	particleSum += _electronDensity * sqrt(eightkT3DivPi / Constant::ELECTRONMASS) * Zsum;

	// (H+) ions
	Zsum = 0;
	for (int z = Zmin; z <= Zmax; z++)
	{
		double ksi = z; // ksi = Ze / q_i = z
		Zsum += stickingCoefficient(a, z, 1) * fZ[z - Zmin] * lambdaTilde(tau, ksi);
	}
	particleSum += _electronDensity * sqrt(eightkT3DivPi / Constant::HMASS_CGS) * Zsum;

	// EA(Zmin) = IP(Zmin-1) because IP(Z) = EA(Z+1)
	double secondTerm = fZ[0] * collisionalChargingRate(a, Zmin, -1, Constant::ELECTRONMASS) *
	                    ionizationPotential(a, Zmin - 1);

	return Constant::PI * a * a * particleSum + secondTerm;
}

void PhotoelectricHeatingRecipe::chargeBalance(double a, const std::vector<double>& wavelengthv,
                                               const std::vector<double>& Qabs,
                                               const std::vector<double>& energyDensity_lambda,
                                               int& resultZmax, int& resultZmin,
                                               std::vector<double>& resultfZ) const
{
	// Express a in angstroms
	double aA = a / Constant::ANG_CM;

	// Shortest wavelength = highest possible energy of a photon
	double hnumax = Constant::PLANCKLIGHT / wavelengthv[0];

	// The maximum charge is one more than the highest charge which still allows ionization by
	// photons of hnumax
	resultZmax = floor(
	                ((hnumax - _workFunction) * Constant::ERG_EV / 14.4 * aA + .5 - .3 / aA) /
	                (1 + .3 / aA));

	// The minimum charge is the most negative charge for which autoionization does not occur
	resultZmin = minimumCharge(a);

	if (resultZmax < resultZmin)
		throw std::range_error("Zmax is smaller than Zmin");
	resultfZ.resize(resultZmax - resultZmin + 1, -1);

	// We will cut off the distribution at some point past the maximum (in either the positive
	// or the negative direction)
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
		double Jpe = emissionRate(a, current, wavelengthv, Qabs, energyDensity_lambda);
		double Jion = collisionalChargingRate(a, current, 1, Constant::HMASS_CGS);
		double Je = collisionalChargingRate(a, current + 1, -1, Constant::ELECTRONMASS);
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
		double Jpe = emissionRate(a, z - 1, wavelengthv, Qabs, energyDensity_lambda);
		double Jion = collisionalChargingRate(a, z - 1, 1, Constant::HMASS_CGS);
		double Je = collisionalChargingRate(a, z, -1, Constant::ELECTRONMASS);
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
		double Jpe = emissionRate(a, z, wavelengthv, Qabs, energyDensity_lambda);
		double Jion = collisionalChargingRate(a, z, 1, Constant::HMASS_CGS);
		double Je = collisionalChargingRate(a, z + 1, -1, Constant::ELECTRONMASS);
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

double PhotoelectricHeatingRecipe::yieldFunctionTest() const
{
	// Parameters
	const int Z = 10;
	vector<double> av = {4e-8, 10e-8, 30e-8, 100e-8, 300e-8};

	// Plot range
	const double hnuMin = 5 / Constant::ERG_EV;
	const double hnuMax = 15 / Constant::ERG_EV;
	const size_t N = 500;

	std::ofstream out("/Users/drvdputt/GasModule/run/yieldTest.dat");
	for (double a : av)
	{
		out << "# a = " << a << '\n';
		const double e2a = Constant::ESQUARE / a;

		// Quantities independent of nu
		double ip_v = ionizationPotential(a, Z);

		double Emin = Z >= 0 ? 0
		                     : -(Z + 1) * e2a /
		                                              (1 +
		                                               std::pow(27. * Constant::ANG_CM / a,
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

double PhotoelectricHeatingRecipe::heatingRateTest(std::string filename) const
{
	readQabs();

	// Wavelength grid
	const vector<double>& frequencyv = Testing::generateGeometricGridv(
	                _nWav, Constant::LIGHT / _maxWav, Constant::LIGHT / _minWav);
	// Input spectrum
	const Array& specificIntensityv = Testing::generateSpecificIntensityv(frequencyv, _Tc, _G0);

	// Convert to wavelength units
	const vector<double>& wavelengthv = Testing::freqToWavGrid(frequencyv);
	Array tempEnergyDensity_lambda =
	                Constant::FPI / Constant::LIGHT *
	                Testing::freqToWavSpecificIntensity(frequencyv, specificIntensityv);
	const vector<double> energyDensity_lambda(begin(tempEnergyDensity_lambda),
	                                          end(tempEnergyDensity_lambda));
	cout << "Made isrf \n";

	//	ofstream spectrumout;
	//	spectrumout.open("/Users/drvdputt/GasModule/run/u_lambda.dat");
	//	for (size_t i = 0; i < wavelengthv.size(); i++)
	//		spectrumout << wavelengthv[i] * Constant::CM_UM << '\t' <<
	//energyDensity_lambda[i] << endl; 	spectrumout.close();

	// Grain sizes for test
	double aMin = 1.5 * Constant::ANG_CM;
	double aMax = 10000 * Constant::ANG_CM;
	const size_t Na = 60;

	// Output file will contain one line for every grain size
	std::ofstream outRate;
	outRate.open(filename);

	//	ofstream outAvgQabs;
	//	outAvgQabs.open("/Users/drvdputt/GasModule/run/avgQabsInterp.txt");

	double aStepFactor = std::pow(aMax / aMin, 1. / Na);
	double a = aMin;

	for (size_t m = 0; m < Na; m++)
	{
		// Absorption spectrum
		vector<double> Qabs = generateQabsv(a, wavelengthv);

		vector<double> QabsTimesUlambda(Qabs.size());
		for (size_t n = 0; n < Qabs.size(); n++)
			QabsTimesUlambda[n] = Qabs[n] * energyDensity_lambda[n];

		double uTimesAvgQabs = NumUtils::integrate<double>(wavelengthv, QabsTimesUlambda);
		double totalAbsorbed = Constant::PI * a * a * Constant::LIGHT * uTimesAvgQabs;

		cout << "Size " << a / Constant::ANG_CM << endl;
		outRate << a / Constant::ANG_CM << '\t'
		        << PhotoelectricHeatingRecipe::heatingRateA(a, wavelengthv, Qabs,
		                                                    energyDensity_lambda) /
		                                totalAbsorbed
		        << '\n';
		//		double u = NumUtils::integrate<double>(wavelengthv,
		//energyDensity_lambda); 		double avgQabs = uTimesAvgQabs / u;
		//		outAvgQabs << a / Constant::ANG_CM << '\t' << avgQabs << endl;
		a *= aStepFactor;
	}
	outRate.close();
	cout << "Wrote " << filename << endl;
	//	outAvgQabs.close();
	cout << "Wrote avgQabsInterp.txt" << endl;

	cout << "Charging parameter = " << _G0 * sqrt(_gasTemperature) / _electronDensity << endl;
	return 0.0;
}

double PhotoelectricHeatingRecipe::chargeBalanceTest() const
{
	readQabs();

// Wavelength grid
#ifdef EXACTGRID
	vector<double> wavelengthv = _filelambdav;
#else
	vector<double> wavelengthv = Testing::generateGeometricGridv(_nWav, _minWav, _maxWav);
#endif

	// Input spectrum
	Array tempIsrf = Testing::generateSpecificIntensityv(wavelengthv, _Tc, _G0);
	vector<double> isrfv(begin(tempIsrf), end(tempIsrf));

	// Grain size
	double a = 200. * Constant::ANG_CM;

	// Absorption spectrum
	vector<double> Qabsv = generateQabsv(a, wavelengthv);

	// Calculate charge distribution
	vector<double> fZv;
	int Zmax, Zmin;
	chargeBalance(a, wavelengthv, Qabsv, isrfv, Zmax, Zmin, fZv);

	std::cout << "Zmax = " << Zmax << " Zmin = " << Zmin << " len fZ = " << fZv.size()
	          << std::endl;

	std::ofstream out;
	out.open("/Users/drvdputt/GasModule/run/fZ.txt");
	out << "# carbon = " << _carbonaceous << endl;
	out << "# a = " << a << endl;
	out << "# Teff = " << _Tc << endl;
	out << "# G0 = " << _G0 << endl;
	out << "# ne = " << _electronDensity << endl;
	out << "# Tgas = " << _gasTemperature << endl;

	for (int z = Zmin; z <= Zmax; z++)
		out << z << '\t' << fZv[z - Zmin] << '\n';
	out.close();

	return 0.;
}
