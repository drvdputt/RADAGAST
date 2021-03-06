#ifndef CORE_LEVELCOEFFICIENTS_HPP
#define CORE_LEVELCOEFFICIENTS_HPP

#include "EigenAliases.hpp"
#include "LevelSolution.hpp"
#include "LineProfile.hpp"
#include "Spectrum.hpp"
#include <array>
#include <memory>

namespace RADAGAST
{
    struct CollisionParameters;

    /** This class deals with level coefficients. A data class (e.g. H2FromFiles) can inherit
        from this class to re-use the level coefficient infrastructure. The most important
        reason for this base class, is a common implemenation of prepareAbsorptionMatrix, which
        integrates over each line to calculate the induced transition rates. T */
    class LevelCoefficients
    {
    public:
        /** An atomic (or molecular) mass needs to be passed, as it will influence the thermal
            broadening of the lines. */
        LevelCoefficients(double mass);
        virtual ~LevelCoefficients() = default;

    protected:
        /** For use by subclass during construction (workaround would be a virtual setup()
            function.) */
        void setConstants(const EVector& ev, const EVector& gv, const EMatrix& avv, const EMatrix& extraAvv);

    public:
        /** Energy of the levels */
        const EVector& ev() const { return _ev; }

        /** Multiplicity of the levels */
        const EVector& gv() const { return _gv; }

        /** Spontaneous radiative transition rates between the levels */
        const EMatrix& avv() const { return _avv; }

        /** Collisional transition rates, calculated by a subclass */
        virtual EMatrix cvv(const CollisionParameters& cp) const = 0;

        /** Spontaneous transition rates which do not produce line emission (really only used for
            two-photon continuum). */
        const EMatrix& extraAvv() const { return _extraAvv; }

        /** Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from
            the Aij) and Pij is the line power, i.e. the radiation field integrated over the
            line profile. The temperature and collision coefficients are needed, because these
            influence the shape of the line profile. The units of Bij and Pij are often
            different in the literature and other codes (it depends on the units used for the
            radiation field), but their product should always have units [s-1]. */
        EMatrix prepareAbsorptionMatrix(const Spectrum& meanIntensity, double T, const EMatrix& Cvv) const;

        /** Ouputs some properties about the different line transitions. The results for the number
            of lines, their frequencies [s-1] and their natural widths (decay rate [s-1] / 4 pi)
            are returned by reference. */
        void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

        /** Calculates the level populations using a simple Boltzman LTE equation. Also serves as
            an example of how to properly set up a LevelSolution object. */
        LevelSolution solveLTE(double density, const CollisionParameters& cp) const;

        LevelSolution solveZero(double T) const;

        /** Return the number of levels in the solution */
        size_t numLv() const { return _numLv; }

        /** The boltzman fractions for the levels, based purely on their energies and the
            temperatures */
        EVector solveBoltzmanEquations(double T) const;

        /** Abstraction of the loop over all lines. Executes thingWithLine for all combinations
            upper > lower that have _Avv(upper, lower) > 0. If the levels are sorted, and all
            downward transitions have line activity, then this function will loop over all elements
            of the lower triangle of the level matrix. */
        void forActiveLinesDo(std::function<void(size_t ini, size_t fin)> thing) const;

        /** Calculates the integrated emission coefficient of a specific line [erg s-1 cm-3 sr-1].
            Multiplying with the line profile [Hz-1] will yield the specific intensity [erg s-1
            cm-3 sr-1 Hz-1]. */
        double lineIntensityFactor(size_t upper, size_t lower, double nu) const;

        /** Computes the integrated opacity of a line [cm-1 Hz]. Multiplying with the line profile
            [Hz-1] will yield the opacity [cm-1] at each frequency. */
        double lineOpacityFactor(size_t upper, size_t lower, double nu, double nl) const;

        /** Return a line profile object that can be used to calculate the (normalized to 1) line
            profile of the "upper-lower" line. The natural line width (the lorenzian contribution)
            is calculated from the total decay rate due to both spontaneous and collisional
            transitions contained in _avv and Cvv, respectively. The thermal line width (the
            gaussian contribution) is calculated from the given temperature and the mass of the
            particle. */
        LineProfile lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const;

    private:
        double _mass;
        size_t _numLv{0};
        EVector _ev;
        EVector _gv;
        EMatrix _avv;
        EMatrix _extraAvv;
    };
}
#endif  // CORE_LEVELCOEFFICIENTS_HPP
