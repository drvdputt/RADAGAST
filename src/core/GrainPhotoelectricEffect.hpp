#ifndef CORE_PHOTOELECTRICHEATING_HPP
#define CORE_PHOTOELECTRICHEATING_HPP

#include "ChargeDistribution.hpp"
#include "Constants.hpp"
#include "GrainInterface.hpp"
#include "Spectrum.hpp"

#include <functional>
#include <vector>

// REMINDER: int i is the size index. These functions used to take double a, but now i should
// refer to an element of sizev, given at construction.

// TODO: Try to fix weird charge balance for some sizes (see code-notes section bugs) when using
// van Hoof's correction (see GasPhysics.pdf)

/** This class provides an implementation of the photoelectric heating recipe described in
    Weingartner \& Draine (2001), hereafter WD01. */
class GrainPhotoelectricCalculator
{
public:
	/** Create a photoelectric effect object. Requires an instance of one of the \c
	GrainType subclasses. The \c SilicateGrain and \c CarbonaceousGrain types are 'native'
	to this recipe, i.e. their implementations are also based on formulae from WD01. With
	some optimism, it should be possible to use this recipe for other types of grains, if
	the right functions and properties (yield, ionization potentials, etc.) are provided. */
	GrainPhotoelectricCalculator(const Array& sizev, double workFunction)
	                : _sizev{sizev}, _workFunction{workFunction}
	{
	}

	double yieldFunctionTest() const;

	/* Makes a plot of the heating efficiency in function of the grain size. Saved as a
	   two-column file in $(pwd)/photoelectric using the filename. */
	void heatingRateTest(double G0, double gasT, double ne) const;

	void chargeBalanceTest(double G0, double gasT, double ne, double np) const;

	/** Gathers the parameters of the gas environment needed to run the photoelectric
	    heating recipe. */
	typedef struct Environment
	{
		/** This constructor creates an environment struct for the photoelectric heating
		    calculation. It takes a frequency grid, a specific intensity of the ambient
		    radiation field for each of those frequencies, a gas temperature, an
		    electron an proton density, and then three lists which specify the particles
		    which participate in the charging of the grains. A list of particle charges,
		    number densities and a list containing the mass of each particle type should
		    be provided. */
		Environment(const Spectrum& specificIntensity, double T, double ne, double np,
		            const std::vector<int>& chargev, const Array& densityv,
		            const Array& massv)
		                : _specificIntensity(specificIntensity), _T(T), _ne(ne),
		                  _np(np), _chargev(chargev), _densityv(densityv),
		                  _massv(massv){};
		// Radiation field
		Spectrum _specificIntensity;
		// The gas temperature, and frequently used densities
		double _T, _ne, _np;
		/* Properties of all the charged particles in the environment (also contains ne
		   and np). */
		std::vector<int> _chargev;
		Array _densityv, _massv;
	} Environment;

	/** Calculates the heating rate per grain for a grain size a. Uses chargeBalance to
	    obtain a charge distribution, and then RateAZ for every charge Z. */
	double heatingRateA(int i, const Environment& env, const Array& Qabsv,
	                    const ChargeDistribution& cd) const;

	/** Uses detailed balance to calculate the charge distribution of a grain a, in and
	    environment env, given the absorption efficiency of that grain in function of the
	    wavelength. */
	ChargeDistribution calculateChargeDistribution(int i, const Environment& env,
	                                               const Array& Qabsv) const;

	/** Recipe from 1991-Baldwin I tried to implement. Returns the cooling per grain surface
	    area [erg s-1 cm-2], so we need to multipy with the surface area per volume unit
	    (n_grain * sigma_grain). The last argument (bool) indicates whether we want to know
	    the cooling of the gas (false), or the heating of the grain (true). In the latter
	    case, it adds a term that takes into account the potential of the grain. */
	double gasGrainCollisionCooling(int i, const Environment& env,
	                                const ChargeDistribution& cd, double Tgrain,
	                                bool addGrainPotential) const;

private:
	/** Implements WD01 equation 24. Calculates the negative charge necessary for a grain to
	    immediately autoionize when an electron is captured. */
	int minimumCharge(int i) const;

	/** Calculates the heating rate by a grain of size a and charge Z, given a
	    wavelength-resolved radiation field and absorption efficiency. */
	double heatingRateAZ(int i, int Z, const Array& frequencyv, const Array& Qabsv,
	                     const Array& specificIntensityv) const;

	/** Calculates the rate at which photoelectrons are emitted from a single grain [s-1],
	    according to equation 25 of WD01. */
	double emissionRate(int i, int Z, const Array& frequencyv, const Array& Qabsv,
	                    const Array& specificIntensityv) const;

	/** The rate [s-1] at which a grain is charged by colliding with other particles. Taken
	    from Draine & Sutin (1987) equations 3.1-3.5. */
	double collisionalChargingRate(int i, double gasT, int Z, int particleCharge,
	                               double particleMass, double particleDensity) const;

	/** The energy removed from the gas by particle sticking to a grain, WD01 equation 42.
	    TODO: figure out how this fits in with the gas-grain collisional energy exchange.
	    I've disabled this for now in Options.hpp. */
	double recombinationCoolingRate(int i, const Environment& env,
	                                const ChargeDistribution& cd) const;

	/** Get the photoelectric and photodetachment thresholds. (Used by both the heating rate
	    integral and number rate integral). */
	void getPET_PDT_Emin(int i, double Z, double& pet, double& pdt, double& Emin) const;

	/** Integrates over the radiation field, counting the number of absorptions per second,
	    per projected grain area. f_hnuDiff can be any function of hnuDiff = (hnu - pet).
	    Multiplying with yield will give the total photoelectric emission rate in electrons
	    s-1 cm-2, while multiplying with the average energy and the yield will give the
	    heating rate in erg s-1 cm-2. */
	double photoelectricIntegrationLoop(
	                const Array& frequencyv, const Array& Qabsv,
	                const Array& specificIntensityv, double pet,
	                const std::function<double(double hnuDiff)>* f_hnuDiff = nullptr) const;

	/** Integration loop which applies equation 20 for the photodetachment cross section. Do
	    not forget to multiply the result with abs(Z)! */
	double
	photodetachmentIntegrationLoop(int Z, const Array& frequencyv,
	                               const Array& specificIntensityv, double pdt,
	                               const double* calcEnergyWithThisEmin = nullptr) const;

	/** Things specific for WD01. If another recipe is ever implemented, then
	    GrainPhotoelectricCalculator should become abstract, and the function and data
	    members indicated here should become virtual/move to a subclass. Each subclass can
	    then have their own caching mechanisms. */
	///@{
	double ionizationPotential(int i, int Z) const;

	double photoelectricYield(int i, int z, double hnuDiff, double Emin) const;

	double autoIonizationThreshold(int i) const;

	double stickingCoefficient(int i, int z, int z_i) const;

	bool _carOrSil;
	///@}

	// Some things will be cached based on the size
	const Array& _sizev;

	// This constant depends on the grain properties. Is copied from the GrainPhotoelectricData
	// object at construction.
	double _workFunction;
};

#endif // CORE_PHOTOELECTRICHEATING_HPP
