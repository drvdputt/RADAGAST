#ifndef _PHOTOELECTRICHEATING_H_
#define _PHOTOELECTRICHEATING_H_

#include "Array.h"
#include "Constants.h"

#include <vector>

// TODO: convert everything to frequency units, or write interface
// TODO: consider van Hoof's correction (see GasPhysics.pdf)
// TODO: create a way to get the excess/deficit of electrons, and maybe charge exchange reaction rates

class GrainPhotoelectricEffect
{
public:
	GrainPhotoelectricEffect(bool carbonaceous)
	                : _carbonaceous(carbonaceous), _workFunction(calcWorkFunction(carbonaceous))
	{
	}

	static double calcWorkFunction(bool carbonaceous)
	{
		return carbonaceous ? 4.4 / Constant::ERG_EV : 8 / Constant::ERG_EV;
	}

private:
	// Grain properties to use for tests (more detailes in generateQabs)
	const bool _carbonaceous;
	const double _workFunction;

	/* Currently, the only public functions are test functions which output some files I can
	   plot. */
public:
	double yieldFunctionTest() const;

	/* Makes a plot of the heating efficiency in function of the grain size. Saved as a
	 two-column file in $(pwd)/photoelectric using the filename. */
	double heatingRateTest(double G0, double gasT, double ne) const;

	double chargeBalanceTest(double G0, double gasT, double ne, double np) const;

private:
	/* Functions to calculate the heating rate according to the recipe by Weingartner and Draine
	   (2001) */

	/* Gathers the parameters of the gas environment which are repeatedly needed. */
	typedef struct Environment
	{
		Environment(const Array& wavv, const Array& u_lambdav, double temperature,
		            double electronDensity, double protonDensity,
		            const std::vector<int>& chargev, const Array& densityv,
		            const Array& massv)
		                : wavelengthv(wavv), energyDensityv(u_lambdav), T(temperature),
		                  ne(electronDensity), np(protonDensity), particleChargev(chargev),
		                  particleDensityv(densityv), particleMassv(massv){};
		// Radiation field
		Array wavelengthv, energyDensityv;
		// The gas temperature, and frequently used densities
		double T, ne, np;
		/* Properties of all the charged particles in the environment (also contains ne and
		   np). */
		std::vector<int> particleChargev;
		Array particleDensityv, particleMassv;
	} Environment;

	/** Calculate the total heating rate by the grains, given a certain environment.  Then, the
	    grain properties need to be specified. First, a list of grain sizes, grainSizev, and a
	    list of number densities for each size, grainDensityv. Then, an absorption efficiency
	    for each grain size and each wavelength should be given. In absQvv, every row represents
	    the wavelength-dependent absorption efficiency for a certain grain size. Therefore, the
	    rows should be indexed in the same way as grainSizev, while the columns should be
	    indexed like env.wavelengthv. Calls heatinRateA for every grain size. */
	double heatingRate(const Environment& env, const std::vector<double>& grainSizev,
	                   const std::vector<double>& grainDensityv,
	                   const std::vector<std::vector<double>>& absQvv) const;

	/** Calculates the heating rate per grain for a grain size a. Uses chargeBalance to obtain a
	    charge distribution, and then calls heatingRateAZ for every charge Z. */
	double heatingRateA(double a, const Environment& env,
	                    const std::vector<double>& Qabs) const;

	/** Uses detailed balance to calculate the charge distribution of a grain a, in and
	    environment env, given the absorption efficiency of that grain in function of the
	    wavelength. */
	void chargeBalance(double a, const Environment& env, const std::vector<double>& Qabs,
	                   int& resultZmax, int& resultZmin, std::vector<double>& resultfZ) const;

	/** Calculates the heating rate by a grain of size a and charge Z, given a
	    wavelength-resolved radiation field and absorption efficiency. */
	double heatingRateAZ(double a, int Z, const Array& wavelengthv,
	                     const std::vector<double>& Qabs,
	                     const Array& energyDensity_lambda) const;

	/** Implements WD01 equations 2, 4 and 5. */
	double ionizationPotential(double a, int Z) const;

	/** Calculates the rate at which photoelectrons are emitted from a single grain [s-1],
	    according to equation 25 of WD01. */
	double emissionRate(double a, int Z, const Array& wavelengthv,
	                    const std::vector<double>& Qabs,
	                    const Array& energyDensity_lambda) const;

	/** Calculates the integral over the energy in WD01 equation 39. */
	double energyIntegral(double Elow, double Ehigh, double Emin, double Emax) const;

	/** Calculates the photoelectric yield according to WD01 equation 12. */
	double yield(double a, int Z, double hnuDiff, double Elow, double Ehigh) const;

	/** Implementes WD01 equations 23, 24. Calculates the negative charge necessary for a grain
	    to immediately autoionize when an electron is captured. */
	int minimumCharge(double a) const;

	/** The sticking coefficient based on equations 1 and 27-30 */
	double stickingCoefficient(double a, int Z, int z_i) const;

	/** The rate [s-1] at which a grain is charged by colliding with other particles. Taken from
	    Draine & Sutin (1987) equations 3.1-3.5. */
	double collisionalChargingRate(double a, double gasT, int Z, int particleCharge,
	                               double particleMass, double particleDensity) const;

	/** Draine & Sutin (1987) equations 3.6-3.10. */
	double lambdaTilde(double tau, double ksi) const;

	/** The energy removed from the gas by particle sticking to a grain, WD01 equation 42. */
	double recombinationCoolingRate(double a, const Environment& env,
	                                const std::vector<double>& fZ, int Zmin) const;

	// Gas properties which worked for photoelectricheatingtest. Keep these as a reference.
	//	const double _hydrogenDensity{2.5e1};
	//	const double _ionizationFraction{3.e-4};
	//	const double _electronDensity{_hydrogenDensity * _ionizationFraction};
	//	double _gasTemperature{1.e3};

	/* The radiation field to use for test will use this blackbody temperature to determine the
	   shape. Its actual strength should be provided using the 'G0' argument of the test
	   function. */
	const double _Tc{3.e4};

	// Wavelength grid to use for tests
	const std::size_t _nWav{200};
	const double _minWav{0.0912 * Constant::UM_CM}; // cutoff at 13.6 eV
	const double _maxWav{1000 * Constant::UM_CM};

	// Data to use for tests
	void readQabs() const;
	std::vector<double> generateQabsv(double a, const std::vector<double>& wavelengthv) const;
};

#endif /* _PHOTOELECTRICHEATING_H_ */
