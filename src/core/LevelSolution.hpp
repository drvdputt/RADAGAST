#ifndef CORE_LEVELSOLUTION_HPP
#define CORE_LEVELSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

namespace RADAGAST
{
    class CollisionParameters;
    class LevelCoefficients;
    class Spectrum;

    /** Storage of non-constant quantities for a level system. Also implements several derived
        quantities. Some data members are just there to keep the memory allocated, and to gather
        all the non-constant level data in one place. Make sure that the arguments passed to the
        setters are compatible are with the LevelCoefficients given at construction. 
        Once all setters have been used correctly, the other functions can be called safely. */
    class LevelSolution
    {
    public:
        /** Pass a pointer to the correct LevelCoefficients object, for access to the
            constants. */
        LevelSolution(const LevelCoefficients* lc) : _levelCoefficients{lc} {};

        /** Set the temperature to T, and the rates and densities to zero, in case the density
            of the relevant species is zero. */
        void setToZero(double T);

        /** Update the stored data, based on the given radiation field and collision parameters,
            using the LevelCoefficients pointer that was given at construction. */
        void updateRates(const Spectrum& meanIntensity, const CollisionParameters& cp);

        /** Same as the above, but with zero radiation field (only collisional) */
        void updateRates(const CollisionParameters& cp);

        /** Return the total transition rate matrix (to be called after updateRates()). This is
            the sum of the spontaneous (Aij), induced (Bij) and collisional (Cij)
            transitions. */
        EMatrix Tvv() const;

        /** Set the level population vector */
        void setNv(const EVector& nv) { _nv = nv; }

        /** Return the stored populations */
        const EVector& nv() const { return _nv; }

        /** Returns true if the level population vector nv is of the wrong dimension, or
            consists entirely of zeros. In that case, it is not suitable for an initial
            guess. */
        bool hasBadNv() const;

        /** Return the fractional populations */
        EVector fv() const;

        /** The spectrum emitted by the line transitions, expressed as the emission coefficient
            j_nu f (erg/cm3/s/Hz). */
        Array emissivityv(const Array& eFrequencyv) const;

        /** The line opacity alpha_nu, equivalent to kappaRho for dust [cm-1]. */
        Array opacityv(const Array& oFrequencyv) const;

        /** Net heating due to (de-)excitation [erg s-1 cm-3] */
        double netHeating() const;

    private:
        const LevelCoefficients* _levelCoefficients;
        double _t{0.};
        EVector _nv;
        EMatrix _bpvv;
        EMatrix _cvv;
    };
}
#endif  // CORE_LEVELSOLUTION_HPP
