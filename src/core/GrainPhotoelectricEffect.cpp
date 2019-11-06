#include "GrainPhotoelectricEffect.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "GrainType.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "WeingartnerDraine2001.hpp"

#include <cmath>
#include <iomanip>

using namespace std;

GrainPhotoelectricEffect::GrainPhotoelectricEffect(const GrainType& grainType)
                : _grainType{grainType}
{
}

int GrainPhotoelectricEffect::minimumCharge(double a) const
{
	double Uait = _grainType.autoIonizationThreshold(a);
	return WD01::minimumCharge(a, Uait);
}

ChargeDistribution
GrainPhotoelectricEffect::calculateChargeDistribution(double a, const Environment& env,
                                                      const Array& Qabsv) const
{
	const Array& frequencyv = env._specificIntensity.frequencyv();
	const Array& specificIntensityv = env._specificIntensity.valuev();

	// Express a in angstroms
	double aA = a / Constant::ANG_CM;

	// Shortest wavelength = highest possible energy of a photon
	double hnumax = Constant::PLANCK * *(end(frequencyv) - 1);

	/* The maximum charge is one more than the highest charge which still allows ionization
	   by photons of hnumax. */
	int resultZmax = floor(
	                ((hnumax - _grainType.workFunction()) * Constant::ERG_EV / 14.4 * aA +
	                 .5 - .3 / aA) /
	                (1 + .3 / aA));

	// The minimum charge is the most negative charge for which autoionization does not
	// occur
	int resultZmin = minimumCharge(a);

	if (resultZmax < resultZmin)
		Error::runtime("Zmax is smaller than Zmin. This can happen when the grains are "
		               "too"
		               "small for the recipe used.");

	// These two functions determine the up and down rates for the detailed balance

	// The rate at which the grain moves out of charge Z, in the positive direction.
	auto chargeUpRate = [&](int Z) -> double {
		double Jtotal{0};
		// Photoelectric effect rate
		Jtotal += emissionRate(a, Z, frequencyv, Qabsv, specificIntensityv);
		/* Collisions with positive particles. Collisions with multiply charged
		   particles are treated as if they can only change the charge by 1. This is
		   technically incorrect, but the charging rate due to positive ions is very
		   small anyway. */
		for (size_t i = 0; i < env._chargev.size(); i++)
			if (env._chargev[i] > 0)
				Jtotal += collisionalChargingRate(a, env._T, Z, env._chargev[i],
				                                  env._massv[i],
				                                  env._densityv[i]);
		return Jtotal;
	};

	// The rate at which the grain moves out of charge Z, in the negative direction.
	auto chargeDownRate = [&](int Z) -> double {
		double Jtotal{0};
		// Collisions with negative particles
		for (size_t i = 0; i < env._chargev.size(); i++)
			if (env._chargev[i] < 0)
				Jtotal += collisionalChargingRate(a, env._T, Z, env._chargev[i],
				                                  env._massv[i],
				                                  env._densityv[i]);
		return Jtotal;
	};

	ChargeDistribution chargeDistribution;
	chargeDistribution.calculateDetailedBalance(chargeUpRate, chargeDownRate, resultZmin,
	                                            resultZmax);

	return chargeDistribution;
}

void GrainPhotoelectricEffect::getPET_PDT_Emin(double a, double Z, double& pet, double& pdt,
                                               double& Emin) const
{
	Emin = WD01::eMin(a, Z);
	double ip_v = _grainType.ionizationPotential(a, Z);
	// WD01 eq 6
	pet = Z >= -1 ? ip_v : ip_v + Emin;
	// Photodetachment, WD01 eq 18 (using EA(Z + 1, a) = IP(Z, a) if Z < 0)
	pdt = ip_v + Emin;
}

double GrainPhotoelectricEffect::photoelectricIntegrationLoop(
                const Array& frequencyv, const Array& Qabsv, const Array& specificIntensityv,
                double pet, const function<double(double hnuDiff)>* f_hnuDiff) const
{
	Array integrandv(frequencyv.size());
	double nu_pet = pet / Constant::PLANCK;

	// yield is zero by definition below photoelectric threshold
	int i = frequencyv.size() - 1;
	while (frequencyv[i] > nu_pet && i > 0)
	{
		double hnu = Constant::PLANCK * frequencyv[i];
		integrandv[i] = Qabsv[i] * specificIntensityv[i] / hnu;
		if (f_hnuDiff)
			integrandv[i] *= (*f_hnuDiff)(hnu - pet);
		i--;
	}
	// <function unit> sr-1 cm-2 s-1
	double integral = TemplatedUtils::integrate<double>(frequencyv, integrandv, max(0, i),
	                                                    frequencyv.size() - 1);
	// <function unit> cm-2 s-1
	return Constant::FPI * integral;
}

double GrainPhotoelectricEffect::photodetachmentIntegrationLoop(
                int Z, const Array& frequencyv, const Array& specificIntensityv, double pdt,
                const double* calcEnergyWithThisEmin) const
{
	Array integrandv(frequencyv.size());
	double nu_pdt = pdt / Constant::PLANCK;

	// no effect below photodetachment threshold
	int i = frequencyv.size() - 1;
	while (frequencyv[i] > nu_pdt && i > 0)
	{
		double hnu = Constant::PLANCK * frequencyv[i];
		double hnuDiff = hnu - pdt;
		// <function unit> / time / angle
		integrandv[i] = WD01::sigmaPDT(Z, hnuDiff) * specificIntensityv[i] / hnu;
		if (calcEnergyWithThisEmin)
			integrandv[i] *= hnuDiff + *calcEnergyWithThisEmin;
		i--;
	}
	// <erg optional> s-1 sr-1
	double integral = TemplatedUtils::integrate<double>(frequencyv, integrandv, max(0, i),
	                                                    frequencyv.size() - 1);
	// <erg optional> s-1
	return Constant::FPI * integral;
}

double GrainPhotoelectricEffect::heatingRateAZ(double a, int Z, const Array& frequencyv,
                                               const Array& Qabsv,
                                               const Array& specificIntensityv) const
{
	const double e2_a = Constant::ESQUARE / a;

	// It's cheaper to calculate these together (and safer, less code duplication), so I
	// made this funtion, and pass these values as arguments when needed.
	double pet, pdt, Emin;
	getPET_PDT_Emin(a, Z, pet, pdt, Emin);

	// WD01 text between eq 10 and 11
	double Elow = Z < 0 ? Emin : -(Z + 1) * e2_a;

	std::function<double(double)> yieldTimesAverageEnergy = [&](double hnuDiff) -> double {
		double Y = _grainType.photoelectricYield(a, Z, hnuDiff, Emin);

		double Ehigh = Z < 0 ? hnuDiff + Emin : hnuDiff;
		// The integral over the electron energy distribution (integral E f(E) dE), over
		// the energy range for which electrons can escape
		double IntE = WD01::energyIntegral(Elow, Ehigh, Emin);
		// Divide by (integral f(E) dE) over the same range, which normalizes the above
		// value --> IntE / y2 gives an average energy
		double y2 = WD01::escapingFraction(Z, Elow, Ehigh);

		return Y * IntE / y2;
	};
	double heatingRatePE =
	                Constant::PI * a * a *
	                photoelectricIntegrationLoop(frequencyv, Qabsv, specificIntensityv, pet,
	                                             &yieldTimesAverageEnergy);

	double heatingRatePD = 0;
	if (Z < 0)
		heatingRatePD = photodetachmentIntegrationLoop(Z, frequencyv,
		                                               specificIntensityv, pdt, &Emin);
	return heatingRatePE + heatingRatePD;
}

double GrainPhotoelectricEffect::heatingRateA(double a, const Environment& env,
                                              const Array& Qabsv,
                                              const ChargeDistribution& cd) const
{
	double totalHeatingForGrainSize = 0;

	for (int Z = cd.zmin(); Z <= cd.zmax(); Z++)
	{
		double fZz = cd.value(Z);
		if (!isfinite(fZz))
			Error::runtime("nan in charge distribution");
		double heatAZ = heatingRateAZ(a, Z, env._specificIntensity.frequencyv(), Qabsv,
		                              env._specificIntensity.valuev());

		/* Fraction of grains in this charge state * heating by a single particle of
		   charge Z. */
		totalHeatingForGrainSize += fZz * heatAZ;
	}

	// The net heating rate (eq 41 without denominator)
	double recCool{0};
	if (Options::grainphotoelectriceffect_recombinationCooling)
		recCool = recombinationCoolingRate(a, env, cd);

	double netHeatingForGrainSize{totalHeatingForGrainSize - recCool};
	return netHeatingForGrainSize;
}

double GrainPhotoelectricEffect::emissionRate(double a, int Z, const Array& frequencyv,
                                              const Array& Qabsv,
                                              const Array& specificIntensityv) const
{
	// Notice that there is quite some duplication compared to heatingRateAZ, but i didn't
	// find it worth the effort to make more abstractions.
	double pet, pdt, Emin;
	getPET_PDT_Emin(a, Z, pet, pdt, Emin);

	function<double(double)> yieldf = [&](double hnuDiff) {
		return _grainType.photoelectricYield(a, Z, hnuDiff, Emin);
	};
	double emissionRatePE = Constant::PI * a * a *
	                        photoelectricIntegrationLoop(frequencyv, Qabsv,
	                                                     specificIntensityv, pet, &yieldf);
	double emissionRatePD = 0;
	if (Z < 0)
		emissionRatePD = photodetachmentIntegrationLoop(
		                Z, frequencyv, specificIntensityv, pdt, nullptr);

	return emissionRatePE + emissionRatePD;
}

double GrainPhotoelectricEffect::collisionalChargingRate(double a, double gasT, int Z,
                                                         int particleCharge,
                                                         double particleMass,
                                                         double particleDensity) const
{
	double kT = Constant::BOLTZMAN * gasT;

	// WD eq 26: akT / q^2 = akT / e^2 / z^2
	double tau = a * kT / Constant::ESQUARE / particleCharge / particleCharge;
	// Ze / q = Z / z
	double ksi = static_cast<double>(Z) / static_cast<double>(particleCharge);

	double Jtilde;
	if (ksi < 0)
	{
		Jtilde = (1. - ksi / tau) * (1. + sqrt(2. / (tau - 2. * ksi)));
	}
	else if (ksi > 0)
	{
		double toSquare = 1. + 1. / sqrt(4. * tau + 3. * ksi);
		Jtilde = toSquare * toSquare * exp(-WD01::thetaKsi(ksi) / tau);
	}
	else
	{
		Jtilde = 1. + sqrt(Constant::PI / 2. / tau);
	}
	return particleDensity * _grainType.stickingCoefficient(a, Z, particleCharge) *
	       sqrt(8. * kT * Constant::PI / particleMass) * a * a * Jtilde;
}

double GrainPhotoelectricEffect::recombinationCoolingRate(double a, const Environment& env,
                                                          const ChargeDistribution& cd) const
{
	// Calculates WD01 equation 42
	double kT = Constant::BOLTZMAN * env._T;
	double eightkT3DivPi = 8 * kT * kT * kT / Constant::PI;

	// For every collision partner, add the contibutions for each possible grain charge.
	double particleSum = 0;
	for (size_t i = 0; i < env._chargev.size(); i++)
	{
		// tau = akT / q^2 (WD01 eq 26)
		int z_i = env._chargev[i];
		double tau = a * kT / z_i / z_i / Constant::ESQUARE;

		double Zsum = cd.sumOverCharge([&](int zGrain) {
			double ksi = zGrain / static_cast<double>(z_i);
			return _grainType.stickingCoefficient(a, zGrain, z_i) *
			       WD01::lambdaTilde(tau, ksi);
		});

		particleSum += env._densityv[i] * sqrt(eightkT3DivPi / env._massv[i]) * Zsum;
	}

	/* The second term of equation 42: autoionization of grains with the most negative
	   charge inhibits the cooling of the gas. */
	// EA(Zmin) = IP(Zmin-1) because IP(Z) = EA(Z+1)
	double secondTerm = 0;
	/* This term is only included when the population of the maximally negative grain charge
	   minimumCharge is significant. If it is not siginicant, then fZ will not cover
	   minimumCharge, (and Zmin > minimumCharge). */
	int zmin = cd.zmin();
	if (zmin == minimumCharge(a))
		secondTerm = cd.value(zmin) *
		             collisionalChargingRate(a, env._T, zmin, -1,
		                                     Constant::ELECTRONMASS, env._ne) *
		             _grainType.ionizationPotential(a, zmin - 1);

	return Constant::PI * a * a * particleSum * kT + secondTerm;
}

double GrainPhotoelectricEffect::gasGrainCollisionCooling(double a, const Environment& env,
                                                          const ChargeDistribution& cd,
                                                          double Tgrain,
                                                          bool addGrainPotential) const
{
	double kT = env._T * Constant::BOLTZMAN;
	double kTgrain = Tgrain * Constant::BOLTZMAN;
	double lambdaG = cd.sumOverCharge([&](int zGrain) {
		double lambdaG_for_this_z = 0;
		double Ug = _grainType.ionizationPotential(a, zGrain);
		double Vg = sqrt(Constant::ESQUARE) * Ug;
		for (size_t i = 0; i < env._massv.size(); i++)
		{
			// Dimensionless
			double ZVg = env._chargev[i] * Vg;
			double psi = ZVg / kT;
			double eta = psi <= 0 ? 1 - psi : exp(-psi);
			double ksi = psi <= 0 ? 1 - psi / 2 : (1 + psi / 2) * exp(-psi);
			double S = _grainType.stickingCoefficient(a, zGrain, env._chargev[i]);
			// cm s-1
			double vbar = sqrt(8 * kT / Constant::PI / env._massv[i]);
			// cm-3 * cm s-1 * erg = cm-2 s-1 erg
			double collisionEnergy = 2 * kT * ksi - 2 * kTgrain * eta;
			if (addGrainPotential)
			{
				// A positively charged particle will slow down before it hits a
				// positively charged grain, and accelerate again when it leaves
				// the grain
				collisionEnergy -= ZVg * eta;
			}
			lambdaG_for_this_z += env._densityv[i] * vbar * S * collisionEnergy;
		}
		return lambdaG_for_this_z;
	});
	// Remember that 1991-Baldwin gives the rate 'per unit projected area'
	return Constant::PI * a * a * lambdaG;
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

		// Quantities independent of nu
		double ip_v = _grainType.ionizationPotential(a, Z);

		double Emin{WD01::eMin(a, Z)};

		// WD01 eq 6
		double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin;

		double hnu = hnuMin;
		const double step = (hnuMax - hnuMin) / N;
		for (size_t n = 0; n < N; n++)
		{
			double hnuDiff = hnu - hnu_pet;
			if (hnuDiff > 0)
				out << hnu * Constant::ERG_EV << '\t'
				    << _grainType.photoelectricYield(a, Z, hnuDiff, Emin)
				    << '\n';
			hnu += step;
		}
		out << '\n';
	}
	out.close();
	return 0.0;
}

namespace
{
// Wavelength grid to use for the tests
void testSpectrum(double G0, Array& frequencyv, Array& specificIntensityv)
{
	const double minWav{0.0912 * Constant::UM_CM}; // cutoff at 13.6 eV
	const double maxWav{1000 * Constant::UM_CM};
	const double Tc{3.e4};
	frequencyv = Testing::generateGeometricGridv(200, Constant::LIGHT / maxWav,
	                                             Constant::LIGHT / minWav);
	specificIntensityv = Testing::generateSpecificIntensityv(frequencyv, Tc, G0);
}
} // namespace

void GrainPhotoelectricEffect::heatingRateTest(double G0, double gasT, double ne) const
{
	Array frequencyv, specificIntensityv;
	testSpectrum(G0, frequencyv, specificIntensityv);

	// Gather environment parameters
	const Spectrum specificIntensity(frequencyv, specificIntensityv);
	const Environment env(specificIntensity, gasT, ne, ne, {-1, 1}, {ne, ne},
	                      {Constant::ELECTRONMASS, Constant::PROTONMASS});

	/* File that writes out the absorption efficiency, averaged using the input radiation
	   field as weights. */
	ofstream avgQabsOf = IOTools::ofstreamFile("photoelectric/avgQabsInterp.txt");

	// Grain sizes for test
	double aMin = 3 * Constant::ANG_CM;
	double aMax = 10000 * Constant::ANG_CM;
	const size_t Na = 90;
	Array sizev = Testing::generateGeometricGridv(Na, aMin, aMax);

	// Output file will contain one line for every grain size
	stringstream efficiencyFnSs;
	efficiencyFnSs << "photoelectric/efficiencyG" << setprecision(4) << scientific << G0
	               << ".dat";
	ofstream efficiencyOf = IOTools::ofstreamFile(efficiencyFnSs.str());

	bool car = true;
	const std::vector<Array> qAbsvv = Testing::qAbsvvForTesting(car, sizev, frequencyv);

	// For every grain size
	for (size_t m = 0; m < Na; m++)
	{
		double a = sizev[m];
		const Array& Qabsv = qAbsvv[m];

		// Integrate over the radiation field
		Array intensityTimesQabsv = Qabsv * specificIntensityv;
		double intensityQabsIntegral = TemplatedUtils::integrate<double>(
		                frequencyv, intensityTimesQabsv);

		// Calculate and write out the charge distribution and heating efficiency
		ChargeDistribution cd = calculateChargeDistribution(a, env, Qabsv);
		stringstream filename;
		filename << "photoelectric/multi-fz/fz_a" << setfill('0') << setw(8)
		         << setprecision(2) << fixed << a / Constant::ANG_CM << ".txt";
		stringstream header;
		header << "# a = " << a << '\n';
		header << "# ne = " << env._ne << '\n';
		header << "# Tgas = " << env._T << '\n';
		cd.plot(filename.str(), header.str());

		double heating = GrainPhotoelectricEffect::heatingRateA(a, env, Qabsv, cd);

		double totalAbsorbed =
		                Constant::PI * a * a * Constant::FPI * intensityQabsIntegral;
		double efficiency = heating / totalAbsorbed;
		if (!isfinite(efficiency))
			cout << "Heating " << heating << " totalabsorbed " << totalAbsorbed
			     << endl;

		efficiencyOf << a / Constant::ANG_CM << '\t' << efficiency << '\n';

		// Calculate and write out the ISRF-averaged absorption efficiency
		double intensityIntegral = TemplatedUtils::integrate<double>(
		                frequencyv, specificIntensityv);
		double avgQabs = intensityQabsIntegral / intensityIntegral;
		avgQabsOf << a / Constant::ANG_CM << '\t' << avgQabs << endl;
	}
	efficiencyOf.close();
	cout << "Wrote " << efficiencyFnSs.str() << endl;
	avgQabsOf.close();
	cout << "Wrote avgQabsInterp.txt" << endl;
	cout << "Charging parameter = " << G0 * sqrt(gasT) / ne << endl;
}

void GrainPhotoelectricEffect::chargeBalanceTest(double G0, double gasT, double ne,
                                                 double np) const
{
	Array frequencyv, specificIntensityv;
	testSpectrum(G0, frequencyv, specificIntensityv);
	Spectrum specificIntensity(frequencyv, specificIntensityv);
	const Environment env(specificIntensity, gasT, ne, np, {-1, 1}, {ne, np},
	                      {Constant::ELECTRONMASS, Constant::PROTONMASS});

	// Grain size
	double a = 200. * Constant::ANG_CM;

	// Qabs for each frequency
	bool car = true;
	Array Qabsv = Testing::qAbsvvForTesting(car, {a}, frequencyv)[0];

	ChargeDistribution cd = calculateChargeDistribution(a, env, Qabsv);

	cout << "Zmax = " << cd.zmax() << " Zmin = " << cd.zmin() << '\n';

	stringstream header;
	header << "# a = " << a << endl;
	header << "# G0 = " << G0 << endl;
	header << "# ne = " << ne << endl;
	header << "# Tgas = " << gasT << endl;
	cd.plot("photoelectric/fZ.txt", header.str());
}
