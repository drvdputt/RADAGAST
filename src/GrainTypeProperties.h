#ifndef GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_
#define GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_

#include "GrainInterface.h"

#include <vector>

/** This namespace groups together some functions that provide grain properties in function of one
    of the grain types. Basically, this is al graintype-specific information, that should actually
    be provided by the grain model implementation of the client code. */
namespace GrainTypeProperties
{
// NEEDED FOR H2 FORMATION //
/** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al.
	   table C.1. */
typedef struct SfcInteractionPar
{
	SfcInteractionPar() = default;
	SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt, double nuH2,
	                  double nuHc, double F);
	/* This boolean is set to false when the default constructor is called, signifiying
		   that the created object does not contain any useful information (meaning that the
		   graintype for which it was contructed is not supported, and that the graintype
		   given in the corresponding static function should be skipped for the H2 formation
		   rate calculation. */
	bool _valid{false};
	const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
} SfcInteractionPar;

/** Return a set of parameters which can be used to calculate the H2 formation rate on the
	    surface of grains of this type. */
SfcInteractionPar sfcInteractionPar(GasModule::GrainType t);

// NEEDED FOR PHOTOELECTRIC HEATING //
/** Indicates wheter a photoelectric heating recipe has been implemented for the given type. */
bool heatingAvailable(GasModule::GrainType t);

/** Work functions for CAR and SIL. */
double workFunction(GasModule::GrainType t);

/** An example implementation which works for CAR and SIL. */
d ouble photoElectricYield(GasModule::GrainType t, double a, int Z, double hnu);

double autoIonizationThreshold(GasModule::GrainType t, double a);

double stickingCoefficient(GasModule::GrainType t, double a, int Z, int z_i);
}
#endif /* GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_ */
