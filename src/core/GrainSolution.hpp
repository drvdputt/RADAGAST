#ifndef CORE_GRAINSOLUTION_HPP
#define CORE_GRAINSOLUTION_HPP

#include "Array.hpp"
#include "ChargeDistribution.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "GrainPopulation.hpp"
#include <vector>

class SpeciesVector;
class Spectrum;

/** This class stores any grain properties that should be adjusted when the densities and
    temperature of the gas change, in relation to a specific GrainPopulation object: a charge
    distribution for every size, an adjusted temperature for every size, and a
    GrainPhotoelectricCalculator which might do some caching. */
class GrainSolution
{
public:
    /** Create a GrainSolution object for the given GrainPopulation. If the given population
        contains a GrainPhotoelectricData instance, then the GrainPhotoelectricCalculator will also
        be initialized. */
    GrainSolution(const GasModule::GrainPopulation* population);

    /** Recalculate the temperature for each size, by finding the temperature for which < blackbody
        * kappa > = < radiation field * kappa > + extra heat (such as collisional, needs to be for
        each size). Some things can probably be cached should this be slow. No idea how to handle
        this if stochastically heated grains would be a thing.

        Calculates several extra contributions to the heating of the grains (collisions
        (Draine and Bertoldi 1996) and heat of H2 formation on the surface (Takahashi
        2001)).

        If this takes up too much time, there is still much room for optimization (e.g.
        reuse some of the values calculated here later, to calculate the heating/cooling of
        the gas.

        Since the heating of the grains due to collisions with gas particles depends on the
        charge distribution of the grains, this should be called after
        recalculateChargeDistributions. */
    void recalculateTemperatures(const GrainPhotoelectricCalculator::Environment& env);

    /** Recalculate the charge distribution for each size using the GrainPhotoelectricCalculator.
        If no calculator is present (population does not have grain photoelectric data), the
        default charge distribution is kept (all grain charges are at 0). */
    void recalculateChargeDistributions(const GrainPhotoelectricCalculator::Environment& env);

    /** Calculate the total energy transfer to the gas due to the thermalization of photoelectrons
        ejected from the grains. Should be called after the charge distributions have been
        calculated (for the same env). */
    double photoelectricGasHeating(const GrainPhotoelectricCalculator::Environment& env) const;

    /** Calculate the total energy transfer from the gas to the grains due to collisions. Should be
        called after the charge distributions have been calculated (for the same env). */
    double collisionalGasCooling(const GrainPhotoelectricCalculator::Environment& env) const;

    /** Calculate the total H2 formation rate per H density unit [s-1], summed over all sizes.
        Should be called after the temperatures grain have been updated. */
    double surfaceH2FormationRateCoeff(double Tgas) const;

private:
    const GasModule::GrainPopulation* _population;

    // Charge distribution for each size
    std::vector<ChargeDistribution> _chargeDistributionv;

    // Adjusted tempererature for each size
    Array _newTemperaturev;

    // Calculator instance which contains caching mechanisms based on the list of sizes of
    // the population.
    std::unique_ptr<GrainPhotoelectricCalculator> _photoelectricCalculator;
};

#endif  // CORE_GRAINSOLUTION_HPP
