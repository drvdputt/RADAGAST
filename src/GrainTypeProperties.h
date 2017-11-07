#ifndef GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_
#define GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_

#include "GrainInterface.h"

#include <vector>

/** This namespace groups together some functions that provide grain properties in function of one
    of the grain types. Basically, this is al graintype-specific information, that should actually
    be provided by the grain model implementation of the client code. */
namespace GrainTypeProperties
{
/** Return a set of parameters which can be used to calculate the H2 formation rate on the
	    surface of grains of this type. */
GasModule::SfcInteractionPar sfcInteractionPar(GasModule::GrainTypeLabel t);

// NEEDED FOR PHOTOELECTRIC HEATING //
/** Indicates wheter a photoelectric heating recipe has been implemented for the given type. */
bool heatingAvailable(GasModule::GrainTypeLabel t);

/** Work functions for CAR and SIL. */
double workFunction(GasModule::GrainTypeLabel t);

/** An example implementation which works for CAR and SIL. */
double photoElectricYield(GasModule::GrainTypeLabel t, double a, int Z, double hnu);

double autoIonizationThreshold(GasModule::GrainTypeLabel t, double a);

double stickingCoefficient(GasModule::GrainTypeLabel t, double a, int Z, int z_i);
}
#endif /* GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_ */
