#ifndef GASMODULE_GIT_SRC_GASGRAININTERACTION_H_
#define GASMODULE_GIT_SRC_GASGRAININTERACTION_H_

namespace GasModule
{
class GrainInterface;
}

namespace GasGrain
{
/** Implementation of the grain surface H2 formation rate recipe decribed by Cazaux & Tielens
    (2002, 2004, 2010) and summarized in Rollig et al 2013 */
double surfaceH2FormationRateCoeff(const GasModule::GrainInterface& g, double Tgas);
} /* namespace GasGrain */
#endif /* GASMODULE_GIT_SRC_GASGRAININTERACTION_H_ */
