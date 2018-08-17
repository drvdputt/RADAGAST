#ifndef GASMODULE_GIT_SRC_GASGRAININTERACTION_H_
#define GASMODULE_GIT_SRC_GASGRAININTERACTION_H_

namespace GasModule
{
class GrainInterface;
}

namespace GasGrain
{
/** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux & Tielens
    (2002, 2004, 2010) and summarized in Rollig et al 2013. This particular implementation
    returns the formation rate without muliplything with nH (atomic hydrogen number density
    [cm-3]). Since nH * rate = [cm-3 s-1], the unit of the returned rate is s-1. */
double surfaceH2FormationRateCoeff(const GasModule::GrainInterface& g, double Tgas);

/** Implementation of a simple gas-grain energy exchange formula. Does not depend on grain
    composition, only works with one dust temperature. From Krumholz 2011 (ApJ, 731, 25). */
double simpleGasGrainCool(double Tdust, double Tgas, double nH, double nH2);

} /* namespace GasGrain */
#endif /* GASMODULE_GIT_SRC_GASGRAININTERACTION_H_ */
