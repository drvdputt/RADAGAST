#ifndef CORE_HMODEL_HPP
#define CORE_HMODEL_HPP

#include "HData.hpp"
#include "LevelSolution.hpp"

namespace GasModule
{
    class HModel
    {
    public:
        /** Pass a pointer to the HData object at construction, so that the constant H data and
            functions can be accessed. A pointer to the radiation field is passed. This way,
            some quantities can be precalculated, using the assumption that the radiation field
            stays constant; with the added benefit that the function singatures become
            smaller. */
        HModel(const HData* hData, const Spectrum* specificIntensity);

        /** Solve the H levels, and store them in this object. */
        void solve(double n, const CollisionParameters& cp);

        /** This function returns the line emission spectrum + the continuum emitted by the 2s-1s
            two-photon process. */
        Array emissivityv(const Array eFrequencyv) const;

        /** This function returns the line opacity. */
        Array opacityv(const Array oFrequencyv) const;

        /** From the level populations, calculate the net heating-cooling balance by
            (de-)excitation */
        double netHeating() const;

        /** The population of the 2s level [cm-3]. This value can be stored to calculate the
            two-photon continuum later, after this object is no longer available, using @c
            TwoPhoton::emissivityv. */
        double n2s() const;

        /** Read acces to the level solution, should one wish to inspect the solution in more
            detail */
        const LevelSolution* levelSolution() const { return &_levelSolution; }

    private:
        /** Returns a vector containing the source terms for the equilibrium equations, i.e. the
            partial recombination rates into each level. Currently, some data from OpenADAS is
            used. [cm-3 s-1] */
        EVector sourcev(const CollisionParameters& cp) const;

        /** Produces the sink term to be used by the equilibrium equations. Currently, some
            hydrogen disappears from the ground state because it's being ionized (ionization
            from higher states is ignored). Any possible sink due to H2 formation is ignored
            since it's much slower than the transitions between which we are trying to calculate
            the balance. [s-1] */
        EVector sinkv() const;

        /** This function calculates the two-photon continuum using Nussbaumer \& Smutz (1984). */
        Array twoPhotonEmissivityv(const Array& eFrequencyv) const;

        const HData* _hData;
        const Spectrum* _specificIntensity;
        double _ionizationRate{0.};
        LevelSolution _levelSolution;
    };
}
#endif  // CORE_HMODEL_HPP
