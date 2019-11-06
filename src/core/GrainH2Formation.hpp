#ifndef CORE_GRAINH2FORMATION_HPP
#define CORE_GRAINH2FORMATION_HPP

#include "Array.hpp"
#include "Constants.hpp"

/** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al. table
    C.1. */
typedef struct SfcInteractionPar
{
	SfcInteractionPar() = default;
	SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt,
	                  double nuH2, double nuHc, double F)
	                : _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt}, _nuH2{nuH2},
	                  _nuHc{nuHc}, _f{F}
	{
	}
	const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
} SfcInteractionPar;

namespace GrainH2FormationData
{
/** Builtin values for this set of parameters for carbonaceous grains */
const SfcInteractionPar carSurface(520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10);

/** Builtin values for this set of parameters for silicate grains */
const SfcInteractionPar silSurface(320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10);

/** Numbers from Takahashi J., Uehara H., 2001, ApJ, 561, 843 for the energy added to a grain
    under H2 formation */
///@{
constexpr double grainHeatingPerH2Formed_sil = 0.4 / Constant::ERG_EV;
constexpr double grainHeatingPerH2Formed_car = 1.72 / Constant::ERG_EV;
///@}
} // namespace GrainH2FormationData

class GrainH2Formation
{
public:
	GrainH2Formation(const SfcInteractionPar& sfcInteractionPar, double heatPerH2)
	                : _sfcInteractionPar{sfcInteractionPar}, _heatPerH2{heatPerH2}
	{
	}

	/** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux &
	    Tielens (2002, 2004, 2010) and summarized in Rollig et al 2013. This particular
	    implementation returns the formation rate without multiplying with nH (atomic
	    hydrogen number density [cm-3]). Since nH * rate = [cm-3 s-1], the unit of the
	    returned rate is s-1. */
	double surfaceH2FormationRateCoeff(const Array& sizev, const Array& temperaturev,
	                                   const Array& densityv, double Tgas) const;

	Array surfaceH2FormationRateCoeffPerSize(const Array& sizev, const Array& temperaturev,
	                                         double Tgas) const;

	Array surfaceH2FormationHeatPerSize(const Array& sizev, const Array& temperaturev,
	                                    double Tgas, double nH) const;

private:
	SfcInteractionPar _sfcInteractionPar;
	double _heatPerH2;
};

#endif // CORE_GRAINH2FORMATION_HPP
