#ifndef _PHOTOELECTRICHEATING_H_
#define _PHOTOELECTRICHEATING_H_

#include "Array.h"
#include "Constants.h"
#include "GrainInterface.h"
#include <functional>
#include <vector>

// TODO: consider van Hoof's correction (see GasPhysics.pdf)
/* TODO: create a way to get the excess/deficit of electrons, and maybe charge exchange reaction
   rates. */

class GrainPhotoelectricEffect
{
public:
	GrainPhotoelectricEffect(const GrainType& grainType);

	// DATA TO USE FOR TESTS

	/* Some static hacks. The first one reads some data ripped from SKIRT into variables
	   declared in the cpp file. Then, a Qabsv for a single size can be generated using
	   generateQabsv. The third function, qAbsvvForTesting does the read step, and then calls
	   the second function for every size of the list. */
	static void readQabs(bool car);
	static Array generateQabsv(double a, const Array& frequencyv);
	static std::vector<Array> qAbsvvForTesting(const Array& av, const Array& frequencyv);

private:
	// TODO remove these
	/* const bool _carbonaceous; */
	/* const double _workFunction; */

	/* Currently, the only public functions are test functions which output some files I can
	   plot. */
public:
	double yieldFunctionTest() const;

	/* Makes a plot of the heating efficiency in function of the grain size. Saved as a
	 two-column file in $(pwd)/photoelectric using the filename. */
	double heatingRateTest(double G0, double gasT, double ne) const;

	double chargeBalanceTest(double G0, double gasT, double ne, double np) const;

private:
	/** Calculates the integral of the partial emission rate (per wavelength interval)
	    multiplied with a certain function of the wavelength f(lambda) of choice. To get the
	    total emission rate (number of electrons per second), one can invoke this function with
	    a lambda function that returns 1. To get the heating rate, f needs to be equal to the
	    average energy of an emitted electron. Since this functions calculates an integral for
	    both photoelectric effect and photodetachtment, two separate functions need to be given.
	    peFunction is the function which will be integrated with the photoelectric yield, while
	    pdFunction will be integrated with the photodissociation cross section. The given
	    functions should take the parameters hnuDiff, Emin and Elow as defined in the body of
	    rateIntegral. */
	double rateIntegral(double a, int Z, double Emin, const Array& frequencyv,
	                    const Array& Qabs, const Array& specificIntensityv,
	                    std::function<double(double hnuDiffpet, double Elow)> peFunction,
	                    std::function<double(double hnuDiffpdt)> pdFunction) const;

public:
	/* Gathers the parameters of the gas environment which are repeatedly needed. */
	typedef struct Environment
	{
		/** This constructor creates an environment struct for the photoelectric heating
		    calculation. It takes a frequency grid, a specific intensity of the ambient
		    radiation field for each of those frequencies, a gas temperature, an electron
		    an proton density, and then three lists which specify the particles which
		    participate in the charging of the grains. A list of particle charges, number
		    densities and a list containing the mass of each particle type should be
		    provided. */
		Environment(const Array& frequencyv, const Array& specificIntensityv, double T,
		            double ne, double np, const std::vector<int>& chargev,
		            const Array& densityv, const Array& massv)
		                : _frequencyv(frequencyv), specificIntensityv(specificIntensityv),
		                  _T(T), _ne(ne), _np(np), _chargev(chargev), _densityv(densityv),
		                  _massv(massv){};
		// Radiation field
		Array _frequencyv, specificIntensityv;
		// The gas temperature, and frequently used densities
		double _T, _ne, _np;
		/* Properties of all the charged particles in the environment (also contains ne and
		   np). */
		std::vector<int> _chargev;
		Array _densityv, _massv;
	} Environment;

	/** Calculate the total heating rate by the grains, given a certain environment.  Then, the
	    grain properties need to be specified. First, a list of grain sizes, grainSizev, and a
	    list of number densities for each size, grainDensityv. Then, an absorption efficiency
	    for each grain size and each wavelength should be given. In absQvv, every row represents
	    the wavelength-dependent absorption efficiency for a certain grain size. Therefore, the
	    rows should be indexed in the same way as grainSizev, while the columns should be
	    indexed like env.wavelengthv. Calls heatinRateA for every grain size. */
	double heatingRate(const Environment&,
	                   const GasModule::GrainInterface::Population& grainPop) const;

private:
	/** Calculates the heating rate per grain for a grain size a. Uses chargeBalance to obtain a
	    charge distribution, and then calls heatingRateAZ for every charge Z. */
	double heatingRateA(double a, const Environment& env, const Array& Qabs) const;

	/** Implements WD01 equation 24. Calculates the negative charge necessary for a grain to
	    immediately autoionize when an electron is captured. */
	int minimumCharge(double a) const;

	/** Uses detailed balance to calculate the charge distribution of a grain a, in and
	    environment env, given the absorption efficiency of that grain in function of the
	    wavelength. */
	void chargeBalance(double a, const Environment& env, const Array& Qabs, int& resultZmax,
	                   int& resultZmin, std::vector<double>& resultfZ) const;

	/** Calculates the heating rate by a grain of size a and charge Z, given a
	    wavelength-resolved radiation field and absorption efficiency. */
	double heatingRateAZ(double a, int Z, const Array& frequencyv, const Array& Qabs,
	                     const Array& specificIntensityv) const;

	/** Calculates the rate at which photoelectrons are emitted from a single grain [s-1],
	    according to equation 25 of WD01. */
	double emissionRate(double a, int Z, const Array& frequencyv, const Array& Qabs,
	                    const Array& specificIntensityv) const;

	/** The rate [s-1] at which a grain is charged by colliding with other particles. Taken from Draine
	    & Sutin (1987) equations 3.1-3.5. */
	double collisionalChargingRate(double a, double gasT, int Z, int particleCharge,
	                               double particleMass, double particleDensity) const;

	/** The energy removed from the gas by particle sticking to a grain, WD01 equation 42. TODO:
	    figure out how this fits in with the gas-grain collisional energy exchange. I've
	    disabled this for now, using the INCLUDERECCOOL macro. */
	double recombinationCoolingRate(double a, const Environment& env,
	                                const std::vector<double>& fZ, int Zmin) const;

	// TODO: THE MEMBERS BELOW SHOULD BE REMOVED EVENTUALLY. THEY ONLY PLAY A ROLE IN SOME OF THE TESTS

	/* The radiation field to use for test will use this blackbody temperature to determine the
	   shape. Its actual strength should be provided using the 'G0' argument of the test
	   function. */
	const double _Tc{3.e4};

	// Wavelength grid to use for tests
	const std::size_t _nWav{200};
	const double _minWav{0.0912 * Constant::UM_CM}; // cutoff at 13.6 eV
	const double _maxWav{1000 * Constant::UM_CM};

	const GrainType& _grainType;
};

#endif /* _PHOTOELECTRICHEATING_H_ */
