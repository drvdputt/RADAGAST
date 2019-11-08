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
	    slow. No idea how to handle this if stochastic heating still needs to work.

	    Calculates several extra contributions to the heating of the grains (collisions
	    (Draine and Bertoldi 1996) and heat of H2 formation on the surface (Takahashi
	    2001)).

	    If this takes up too much time, there is still much room for optimization (e.g.
	    reuse some of the values calculated here later, to calculate the heating/cooling of
	    the gas. */
	void recalculateTemperatures(const GrainPhotoelectricCalculator::Environment& env);

	void
	recalculateChargeDistributions(const GrainPhotoelectricCalculator::Environment& env);
	double
	photoelectricGasHeating(const GrainPhotoelectricCalculator::Environment& env) const;
	double
	collisionalGasCooling(const GrainPhotoelectricCalculator::Environment& env) const;

	/** H2 formation rate which depends on the grain temperature solution */
	double surfaceH2FormationRateCoeff(double Tgas) const;

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
