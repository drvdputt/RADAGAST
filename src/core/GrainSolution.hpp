#ifndef CORE_GRAINSOLUTION_HPP
#define CORE_GRAINSOLUTION_HPP

#include "Array.hpp"
#include "ChargeDistribution.hpp"
#include "GasDiagnostics.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "GrainPopulation.hpp"
#include "Testing.hpp"
#include <vector>

namespace RADAGAST
{
    class SpeciesVector;
    class Spectrum;

    /** This class stores any grain properties that should be adjusted when the densities and
        temperature of the gas change, in relation to a specific GrainPopulation object: a charge
        distribution for every size, an adjusted temperature for every size, and a
        GrainPhotoelectricCalculator which might do some caching. */
    class GrainSolution
    {
        friend void Testing::writeGrains(const std::string& outputPath, const std::vector<GrainSolution>& grs,
                                         bool bulkCar);

    public:
        /** Create a GrainSolution object for the given GrainPopulation. If the given population
            contains a GrainPhotoelectricData instance, then the GrainPhotoelectricCalculator
            will also be initialized. The radiation field is assumed to be constant during the
            lifetime of this object, so pass a pointer to it here. */
        GrainSolution(const RADAGAST::GrainPopulation* population, const Spectrum* meanIntensity);

        /** Recalculate the grain temperatures and grain charges using radiation field given at
            construction, for new temperature and densities. */
        void recalculate(double T, const SpeciesVector& sv);

    private:
        /** Recalculate the temperature for each size, by finding the temperature for which <
            blackbody kappa > = < radiation field * kappa > + extra heat (such as collisional,
            needs to be for each size). Some things can probably be cached should this be slow. No
            idea how to handle this if stochastically heated grains would be a thing.

            Calculates several extra contributions to the heating of the grains (collisions (Draine
            and Bertoldi 1996) and heat of H2 formation on the surface (Takahashi 2001)).

            If this takes up too much time, there is still much room for optimization (e.g. reuse
            some of the values calculated here later, to calculate the heating/cooling of the gas.

            Since the heating of the grains due to collisions with gas particles depends on the
            charge distribution of the grains, this should be called after
            recalculateChargeDistributions. */
        void recalculateTemperatures();

        /** Recalculate the charge distribution for each size using the
            GrainPhotoelectricCalculator. If no calculator is present (population does not have
            grain photoelectric data), the default charge distribution is kept (all grain charges
            are at 0). */
        void recalculateChargeDistributions();

    public:
        /** Calculate the total energy transfer to the gas due to the thermalization of
            photoelectrons ejected from the grains. Should be called after @c recalculate(). */
        double photoelectricGasHeating();

        /** Calculate the total energy transfer from the gas to the grains due to collisions.
            Should be called after @c recalculate(). */
        double collisionalGasCooling() const;

        /** Calculate the total H2 formation rate per H density unit [s-1], summed over all sizes.
            Should be called after @c recalculate(). */
        double surfaceH2FormationRateCoeff(double Tgas) const;

        /** Return the average charge for grains with size index m of this population. Assumes
            that recalculate() has already been called. */
        double averageCharge(int m) const;

        /** Add some info about the grains to the diagnostics. Choose a unique prefix for the
            diagnostic names to prevent output conflicts, because multiple populations are
            possible. */
        void extraDiagnostics(GasDiagnostics& gd, const std::string& prefix) const;

    private:
        // pointers to constant data
        const RADAGAST::GrainPopulation* _population;
        const Spectrum* _meanIntensity;
        int _numSizes;

        // Calculator instance which contains caching mechanisms based on the list of sizes of the
        // population.
        std::unique_ptr<GrainPhotoelectricCalculator> _photoelectricCalculator;

        // Workspace for the calculator
        GrainPhotoelectricCalculator::Locals _photoelectricLocals;

        // Charge distribution for each size
        std::vector<ChargeDistribution> _chargeDistributionv;

        // Adjusted temperature for each size
        Array _newTemperaturev;

        // Heating (of grain) due to H2 formation for each size
        Array _h2Heatv;

        // Cache photoelectric heating (or gas) during recalculate(), so it can be reused for
        // photoelectricGasHeating. This is the heating for a single grain of a certain size.
        Array _cachedPEHeatv;
    };
}
#endif  // CORE_GRAINSOLUTION_HPP
