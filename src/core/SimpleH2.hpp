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
        SimpleH2(const LookupTable* lteCool);

        /** Calculate the H2 to H2* ratio, and all quantities that depend on the radiation
            field. Needs to be called before any of the other functions. */
        void solve(double n, const CollisionParameters& cp, const Spectrum& specificIntensity,
                   double h2form = 0) override;
        double dissociationRate(const Spectrum& specificIntensity) const override;
        double dissociationHeating(const Spectrum& specificIntensity) const override;
        double netHeating() const override;
        double orthoPara() const override;
        Array emissivityv(const Array& eFrequencyv) const override;
        Array opacityv(const Array& oFrequencyv) const override;

    private:
        double gloverAbel08Cooling(const CollisionParameters& cp) const;

        // Question: are we going to have H2 and H2* separately in the chemical network? In that
        // case, the big H2 model also needs a way to convert the level populations in to H2 and
        // H2* densities (shouldn't be hard, just use an energy threshold somewhere).
        double _nH2{0};             // total H2 density
        double _nH2g{0};            // H2 density in pseudo ground state
        double _nH2s{0};            // H2* density (vib-rot excited pseudo level)
        double _ortho{.75};         // ortho fraction
        double _para{1. - _ortho};  // para fraction
        double _g{0};               // radiation field in habing units
        double _solomonDissoc{0};   // eq. A8
        double _directDissoc{0};    // eq. A12
        double _gamma3{0};          // Heating due to solomon dissociation, eq. A9
        // Collisional de-excitation heating, eq. A13, due to UV pumping (first to another E
        // level, then back to X, with higher J or v)
        double _gamma4{0};
        // Cooling due to upward excitation by H, H2, e-, p+.
        double _collisionalExcitationCooling{0};
        // Cooling in the high density limit (LTE, precalculated).
        const LookupTable* _lteCool;
    };
}
#endif  // CORE_SIMPLEH2_HPP
