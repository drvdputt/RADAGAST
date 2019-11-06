#ifndef CORE_GRAINSOLUTION_HPP
#define CORE_GRAINSOLUTION_HPP

#include "Array.hpp"
#include "ChargeDistribution.hpp"

#include <vector>

class GrainPhotoelectricCalculator;
class GrainPopulation;
class SpeciesVector;
class Spectrum;

class GrainSolution
{
public:
	GrainSolution(const GrainPopulation& population);

	/** Recalculate the temperature for each size, by finding the temperature for which <
	    blackbody * kappa > = < radiation field * kappa > + extra heat (such as collisional,
	    needs to be given for each size). Some things can probably be cached should this be
	    slow. No idea how to handle this if stochastic heating still needs to work. */
	void recalculateTemperatures(const Spectrum& specificIntensity, Array otherGrainHeat);

	void recalculateChargeDistributions(const Spectrum& specificIntensity,
	                                    const SpeciesVector& speciesNv);
	double photoelectricGasHeating(const SpeciesVector& speciesNv);
	double collisionalGasCooling(const SpeciesVector& speciesNv);
	double collisionalGrainHeating(const SpeciesVector& speciesNv);

private:
	const GrainPopulation* _population;

	// Charge distribution for each size
	std::vector<ChargeDistribution> _chargeDistributionv;

	// Adjusted tempererature for each size
	Array _newTemperaturev;

	// Calculator instance which contains caching mechanisms based on the list of sizes of
	// the population.
	std::unique_ptr<GrainPhotoelectricCalculator> _photoelectricCalculator;
};

#endif // CORE_GRAINSOLUTION_HPP
