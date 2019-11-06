#ifndef CORE_GRAINSOLUTION_HPP
#define CORE_GRAINSOLUTION_HPP

#include "Array.hpp"
#include "ChargeDistribution.hpp"
#include "GrainPhotoelectricEffect.hpp"

#include <vector>

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

	void recalculateChargeDistributions(const GrainPhotoelectricCalculator::Environment& env);
	double photoelectricGasHeating(const GrainPhotoelectricCalculator::Environment& env) const;
	double collisionalGasCooling(const GrainPhotoelectricCalculator::Environment& env) const;

private:
	const GrainPopulation* _population;

	// Charge distribution for each size
	std::vector<ChargeDistribution> _chargeDistributionv;

	// Adjusted tempererature for each size
	Array _newTemperaturev;

	// Calculator instance which contains caching mechanisms based on the list of sizes of
	// the population.
	std::unique_ptr<GrainPhotoelectricCalculator> _photoelectricCalculator{nullptr};
};

#endif // CORE_GRAINSOLUTION_HPP
