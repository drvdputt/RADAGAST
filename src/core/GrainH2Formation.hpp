#ifndef CORE_GRAINH2FORMATION_HPP
#define CORE_GRAINH2FORMATION_HPP

#include "Array.hpp"
#include "Constants.hpp"

namespace RADAGAST
{
    /** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al.
        table C.1. */
    typedef struct SfcInteractionPar
    {
        SfcInteractionPar() = default;
        SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt, double nuH2, double nuHc,
                          double F)
            : _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt}, _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
        {}
        const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
    } SfcInteractionPar;

    namespace GrainH2FormationData
    {
        /** Builtin values for this set of parameters for carbonaceous grains */
        const SfcInteractionPar carSurface(520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10);

        /** Builtin values for this set of parameters for silicate grains */
        const SfcInteractionPar silSurface(320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10);

        /** Numbers from Takahashi J., Uehara H., 2001, ApJ, 561, 843 for the energy added to a
            grain under H2 formation */
        ///@{
        constexpr double grainHeatingPerH2Formed_sil = 0.4 * Constant::EV;
        constexpr double grainHeatingPerH2Formed_car = 1.72 * Constant::EV;
        ///@}
    }  // namespace GrainH2FormationData

    /** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux & Tielens
        (2002, 2004, 2010; check the errata!) and summarized in Rollig et al 2013. */
    class GrainH2Formation
    {
    public:
        /** Create a new GrainH2Formation object, storing a set of parameters to describe the
            formation rate, and an average amount of energy deposited per H2 formed [erg]. */
        GrainH2Formation(const SfcInteractionPar& sfcInteractionPar, double heatPerH2)
            : _sfcInteractionPar{sfcInteractionPar}, _heatPerH2{heatPerH2}
        {}

        /** This particular implementation returns the formation rate without multiplying with nH
            (atomic hydrogen number density [cm-3]). Since nH * rate = [cm-3 s-1], the unit of the
            returned rate is s-1. Three arrays need to be passed: one with the grain sizes, one
            with the grain temperatures, and one with the grain densities. */
        double surfaceH2FormationRateCoeff(const Array& sizev, const Array& temperaturev, const Array& densityv,
                                           double Tgas) const;

        /** Get the H2 formation coefficients as an array, containing one element per grain size.
            [s-1] */
        Array surfaceH2FormationRateCoeffPerSize(const Array& sizev, const Array& temperaturev, double Tgas) const;

        /** Get the heat deposited into the grains, due to H2 formation on their surfaces. This is
            calculated individually per grain size. [erg s-1] */
        Array surfaceH2FormationHeatPerSize(const Array& sizev, const Array& temperaturev, double Tgas,
                                            double nH) const;

    private:
        SfcInteractionPar _sfcInteractionPar;
        double _heatPerH2;
    };
}
#endif  // CORE_GRAINH2FORMATION_HPP
