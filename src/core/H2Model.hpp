#ifndef CORE_H2MODEL_HPP
#define CORE_H2MODEL_HPP

#include "Array.hpp"
#include "GasDiagnostics.hpp"
#include "LevelSolution.hpp"

namespace RADAGAST
{
    struct CollisionParameters;
    class Spectrum;

    /** Abstract class for the non-constant H2 properties, and the functions that depend on
        them. We envision 2 subclasses: a large model which includes a full calculation of the
        level populations, and dissociation / heating / etc. derived from these populations; and
        a smaller model which implements more approximate recipes. */
    class H2Model
    {
    public:
        virtual ~H2Model() = default;

        /** Solve the state of H2, given a constant set of environmental parameters. Depending
            on the sublass, the calculation can be trivial, or extremely heavy (solving the
            levels). The H2 formation tempo [cm-3 s-1] (from grain surfaces only) is need to
            normalize the source term, which describes the formation pumping. It is assumed that
            the radiation field remains constant during the lifetime of the H2 model, hence this
            function does not take it as an argument. Subclasses can calculate the necessary
            quantities derived from the radiation field at construction or in a specialized
            setter. */
        virtual void solve(double n, const CollisionParameters& cp, double h2form = 0) = 0;

        /** The dissociation rate, both by direct photodissociation and the indirect Solomon
            process. This rate can be used to calculate the H2 abundance using a chemical network.
            [s-1] */
        virtual double dissociationRate() const = 0;

        /** Currently takes into account the leftover kinetic energy after direct dissociation
            (effect similar to H ionization), and after Solomon dissociation. */
        virtual double dissociationHeating() const = 0;

        /** Calculate the net heating-cooling balance by (de-)excitation. Typically, the gas is
            heated by UV pumping of H2 (X + UV photon -> B or C -> X' -> X + heat), but cooled by
            the upward collisions (X + heat -> X') */
        virtual double netHeating() const = 0;

        /** Ortho-para ratio (ortho / total) */
        virtual double orthoPara() const = 0;

        /** H2 emission (if any) */
        virtual Array emissivityv(const Array& eFrequencyv) const = 0;

        /** H2 opacity, both by lines and the continuum dissociation cross section (if any) */
        virtual Array opacityv(const Array& oFrequencyv) const = 0;

        /** For debug / diagnostics purposes, a pointer to the level solution object can be given,
            if one exists in the subclass implementation */
        virtual bool hasLevels() const { return false; }
        virtual const LevelSolution* levelSolution() const { return nullptr; }

        /** Add extra diagnostics related to H2 to the given GasDiagnostics object. Does nothing
            by default, but can be overridden by the big or small H2 model to add different
            contributions (preferably things that aren't already covered by
            GasSolution::fillDiagnostics). */
        virtual void extraDiagnostics(GasDiagnostics&) const {}
    };
}
#endif  // CORE_H2MODEL_HPP
