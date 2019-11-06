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
	                  double nuH2, double nuHc, double F);
	/* This boolean is set to false when the default constructor is called, signifiying that
	   the created object does not contain any useful information (meaning that the
	   graintype for which it was contructed is not supported, and that the graintype given
	   in the corresponding static function should be skipped for the H2 formation rate
	   calculation. */
	bool _valid{false};
	const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
} SfcInteractionPar;

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

class GrainH2Formation
{
public:
	GrainH2Formation(const SfcInteractionPar& sfcInteractionPar, double heatPerH2);

	/** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux &
	    Tielens (2002, 2004, 2010) and summarized in Rollig et al 2013. This particular
	    implementation returns the formation rate without multiplying with nH (atomic
	    hydrogen number density [cm-3]). Since nH * rate = [cm-3 s-1], the unit of the
	    returned rate is s-1. */
	double surfaceH2FormationRateCoeff(const Array& sizev, double Tgas);

	Array surfaceH2FormationRateCoeffPerSize(const Array& sizev, double Tgas);

	Array surfaceH2FormationHeatPerSize(const Array& sizev, double Tgas, double nH);

private:
	SfcInteractionPar _sfcInteractionPar;
	double _heatPerH2;
};

#endif // CORE_GRAINH2FORMATION_HPP
