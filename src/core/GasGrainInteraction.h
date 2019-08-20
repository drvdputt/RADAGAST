#ifndef GASMODULE_GIT_SRC_GASGRAININTERACTION_H_
#define GASMODULE_GIT_SRC_GASGRAININTERACTION_H_

#include "Array.h"
#include "Constants.h"
#include "GrainInterface.h"

/** Originally, the idea was to put all direct gas-grain interaction formulae here. But since
    the formula for the gas-grain collisions depends on the grain charges, I have moved that one
    to the grain photoelectric effect namespace. For now, I will only put grain surface related
    things in here. */
namespace GasGrain
{
/** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux & Tielens
    (2002, 2004, 2010) and summarized in Rollig et al 2013. This particular implementation
    returns the formation rate without muliplything with nH (atomic hydrogen number density
    [cm-3]). Since nH * rate = [cm-3 s-1], the unit of the returned rate is s-1. */
double surfaceH2FormationRateCoeff(const GasModule::GrainInterface& g, double Tgas);

Array surfaceH2FormationRateCoeffPerSize(const GasModule::GrainInterface::Population& pop,
                                         double Tgas);

Array surfaceH2FormationHeatPerSize(const GasModule::GrainInterface::Population& pop,
                                    double Tgas, double nH);

/** Numbers from Takahashi J., Uehara H., 2001, ApJ, 561, 843 for the energy added to a grain
    under H2 formation */
constexpr double grainHeatingPerH2Formed_sil = 0.4 / Constant::ERG_EV;
constexpr double grainHeatingPerH2Formed_car = 1.72 / Constant::ERG_EV;

const GasModule::SfcInteractionPar carSurface(520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10);
const GasModule::SfcInteractionPar silSurface(320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10);

} /* namespace GasGrain */
#endif /* GASMODULE_GIT_SRC_GASGRAININTERACTION_H_ */
