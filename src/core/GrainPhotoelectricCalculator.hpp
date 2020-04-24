#ifndef CORE_GRAINPHOTOELECTRICCALCULATOR_HPP
#define CORE_GRAINPHOTOELECTRICCALCULATOR_HPP

#include "ChargeDistribution.hpp"
#include "Constants.hpp"
#include "GrainInterface.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include <array>
#include <functional>
#include <vector>

namespace GasModule
{
    /** This class provides an implementation of the photoelectric heating recipe described in
        Weingartner \& Draine (2001), hereafter WD01. There is some caching of grain properties
        to speed up the calculation. Because in a simulation, the grain properties could be
        different per cell, a new instance of this object is created every time the gas is
        updated. (Reminder: GasSolution has multiple GrainSolutions, and each GrainSolution has
        a PhotoelectricCalculator.) It should be OK to have a changing state in this class,
        which could be very useful for more advanced caching. Many functions are not const for
        that reason.

        On the other hand, I could force the user to make sure that their GrainInterface is
        thread local, and then move this class, keeping the caching mechanisms, into
        GrainInterface. This would make it possible to keep the cache between cell updates. But
        let's try keeping it as it is for now, and seeing if I can speed things up by caching
        some of the more expensive integrands (and hoping the allocations to do the caching
        don't cause to much trouble).

        With some optimism, it should be possible to use this recipe for other types of grains,
        if the right functions and properties (yield, ionization potentials, etc.) are made
        abstract, and implemented in specific subclasses per grain type. */
    class GrainPhotoelectricCalculator
    {
    public:
        /** Currently, this constructor also takes arguments related to the WD01 grain types. This
            should change if the GrainPhotoelectricData and GrainPhotoelectricCalculator classes
            are ever made abstract, to acommodate for other grain types. */
        GrainPhotoelectricCalculator(const Array* sizev, const std::vector<Array>* qAbsvv, double workFunction,
                                     bool carOrSil, const Spectrum* meanIntensity);

        /** Gathers the parameters that are different depending on the iteration, or the thread
            using this calculator. */
        class Locals
        {
            friend class GrainPhotoelectricCalculator;

        public:
            Locals() = default;

            /** This constructor creates an environment struct for the photoelectric heating
                calculation. It takes the gas temperature, and the species densities. From the
                latter, the most important densities are extracted (ne, np, nH, nH2), and those
                densities are placed in an array at the same indices as the charges and masses,
                for easy iteration over the relevant species. */
            Locals(double T, const SpeciesVector& sv);

        private:
            double _T{0.};
            static const int _ie, _ip, _iH, _iH2;
            static const std::array<int, 4> _chargev;
            static const std::array<double, 4> _massv;
            std::array<double, 4> _densityv;
        };

        /** Calculates the heating rate per grain for a grain size a. Uses chargeBalance to
            obtain a charge distribution, and then RateAZ for every charge Z. */
        double heatingRateA(int i, const ChargeDistribution& cd);

        /** Uses detailed balance to calculate the charge distribution of a grain a, in and
            environment env, given the absorption efficiency of that grain in function of the
            wavelength. */
        void calculateChargeDistribution(int i, Locals& env, ChargeDistribution& cd);

        /** The cooling due to collisions with a single grain of the given size [erg s-1]. This
            function can also be used to calculate the extra heat that goes into the grain because
            of this process. To do this, the last argument can be set to @c true, which adds two
            extra terms (see recipe from 1991-Baldwin), which take into account the potential of
            the grain and the energy of recombinations happening on its surface. */
        double gasGrainCollisionCooling(int i, const Locals& env, const ChargeDistribution& cd, double Tgrain,
                                        bool forGrain) const;

        /** The energy removed from the gas due to charged particles recombining with a grain, WD01
            equation 42. This is disabled in Options.hpp, because I think this conflicts (as in,
            double counts something) with the gas-grain collision recipe. */
        double recombinationCoolingRate(int i, const Locals& env, const ChargeDistribution& cd) const;

    private:
        /** Implements WD01 equation 24. Calculates the negative charge necessary for a grain to
            immediately autoionize when an electron is captured. */
        int minimumCharge(int i) const;

        /** Calculates the heating rate by a grain of size a and charge Z, given its absorption
            efficiency. */
        double heatingRateAZ(int i, int Z);

        /** Calculates the rate at which photoelectrons are emitted from a single grain [s-1],
            according to equation 25 of WD01. */
        double emissionRate(int i, int Z);

        /** The rate [s-1] at which a grain is charged by colliding with other particles. Taken
            from Draine & Sutin (1987) equations 3.1-3.5. */
        double collisionalChargingRate(int i, double gasT, int Z, int particleCharge, double particleMass,
                                       double particleDensity) const;

        /** Get the photoelectric and photodetachment thresholds. (Used by both the heating rate
            integral and number rate integral). */
        void getPET_PDT_Emin(int i, int Z, double& pet, double& pdt, double& Emin) const;

        /** Integrates over the radiation field, counting the number of absorptions per second,
            per projected grain area, times a given tabulated function fNu. The integration over
            frequency is only carried about above the photoelectric threshold. The values in fNu
            below the index of the photoelectric threshold do not matter (index referes to
            position in env._meanIntensity.frequencyv()). When fNu is the yield for each
            frequency, this function will give the total photoelectric emission rate in
            electrons s-1 cm-2, while passing the average energy times yield will give the
            heating rate in erg s-1 cm-2. fNu needs to be the same size as
            env._meanIntensity. */
        double photoelectricIntegrationLoop(int i, double nuPET, const Array& fNu);

        /** Integration loop which applies equation 20 for the photodetachment cross section. Do
            not forget to multiply the result with abs(Z)! */
        double photodetachmentIntegrationLoop(int Z, double pdt, const double* calcEnergyWithThisEmin = nullptr);

        /** Things specific for WD01. If another recipe is ever implemented, then
            GrainPhotoelectricCalculator should become abstract, and the function and data members
            indicated here should become virtual/move to a subclass. Each subclass can then have
            their own caching mechanisms. */
        ///@{
        double ionizationPotential(int i, int z) const;

        double photoelectricYield(int i, int z, double hnuDiff, double Emin) const;

        double autoIonizationThreshold(int i) const;

        double stickingCoefficient(int i, int z, int z_i) const;
        ///@}

        // constant data
        const Spectrum* _meanIntensity;
        double _workFunction;
        bool _carOrSil;
        const Array* _sizev;
        const std::vector<Array>* _qAbsvv;

        // constant after construction (mostly precalculated data as function of grain size)
        Array _y1Cache;  // output of WD01::y1 (contains expm1)
        // sticking coefficients for electrons colliding with positive and negative grains (both
        // contain expm1)
        Array _eStickPositiveCache;
        Array _eStickNegativeCache;

        // cached data (changes during the calculation)
        Array _integrationWorkspace;  // some memory to keep allocated between integrations

        // Cache stuff for pairs of i (size index), z (charge). Use map for now, but in
        // principle, a contiguous array could be used, since we have a fixed upper limit for
        // the amount of charges anyway.

        // when the emission rate has been calculated, we can reuse the yield (as function of
        // frequency) for the heating rate calculation. Might get rid of this later, if it's not
        // worth the effort.
        std::map<std::array<int, 2>, Array> _yieldCache;

        // as long as the grain properties and the radiation field don't change, these only
        // depend on the charge and the size. By caching these, new values will only have to be
        // calculated if the range of one of the charge distributions changes (which can happen
        // when gas temperature changes)
        std::map<std::array<int, 2>, double> _emissionRateCache;
        std::map<std::array<int, 2>, double> _heatingRateCache;

    public:
        // some utilities for testing
    };
}
#endif  // CORE_GRAINPHOTOELECTRICCALCULATOR_HPP
