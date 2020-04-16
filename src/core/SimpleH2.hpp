#ifndef CORE_SIMPLEH2_HPP
#define CORE_SIMPLEH2_HPP

#include "H2Model.hpp"

namespace GasModule
{
    class LookupTable;

    /** This class implement a simple model for H2, based mostly on the H2 model of Tielens and
        Hollenbach (1985). For the cooling by H2 lines, cooling curves from Glover and Abel (2008)
        are used instead. Since this is a light weight object, pointers to some externally loaded
        data need to be provided at construction. */
    class SimpleH2 : public H2Model
    {
    public:
        /** Create a new workspace for the simple H2 model. A pointer to the LTE H2 line cooling
            table needs to be given (see SpeciesModelManager). */
        SimpleH2(const LookupTable* lteCool, const Spectrum* specificIntensity);

        /** Calculate the H2 to H2* ratio, and all other quantities that depend on more than
            just the radiation field. Should be called before any of the functions below. */
        void solve(double n, const CollisionParameters& cp, double h2form = 0) override;
        double dissociationRate() const override;
        double dissociationHeating() const override;
        double netHeating() const override;
        double orthoPara() const override;
        Array emissivityv(const Array& eFrequencyv) const override;
        Array opacityv(const Array& oFrequencyv) const override;

    private:
        double gloverAbel08Cooling(const CollisionParameters& cp) const;

        // set at construction
        double _g{0};  // radiation field in habing units

        // set during solve()
        double _nH2{0};                           // total H2 density
        double _nH2g{0};                          // H2 density in pseudo ground state
        double _nH2s{0};                          // H2* density (vib-rot excited pseudo level)
        double _ortho{.75};                       // ortho fraction
        double _para{1. - _ortho};                // para fraction
        double _solomonDissoc{0};                 // eq. A8
        double _directDissoc{0};                  // eq. A12
        double _gamma3{0};                        // heating due to solomon dissociation, eq. A9
        double _gamma4{0};                        // UV pump + collisional de-excitation heating, eq. A13
        double _collisionalExcitationCooling{0};  // cooling due to upward excitation by H, H2, e-, p+.

        // Cooling in the high density limit (LTE, precalculated).
        const LookupTable* _lteCool;

        // Note: if H2 and H2* are treated separately in the chemical network, we will need a
        // way to calculate this ratio for the big H2 model too
    };
}
#endif  // CORE_SIMPLEH2_HPP
