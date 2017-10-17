#ifndef GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_
#define GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_

#include "Array.h"
#include "GrainInterface.h"
#include "GrainPhotoelectricEffect.h"

#include <vector>

/** This class groups together some functions that provide grain properties in function of one of
    the grain types. Most of the functions in here will be static, so maybe this should just be a
    namespace. */
class GrainTypeProperties
{
public:
	GrainTypeProperties();

	/** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al.
	   table C.1. */
	typedef struct SfcInteractionPar
	{
		SfcInteractionPar() = default;
		SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt,
		                  double nuH2, double nuHc, double F);
		/* This boolean is set to false when the default constructor is called, signifiying
		   that no number are available, and that the graintype given in the corresponding
		   static function should be skipped for the H2 formation rate calculation. */
		bool _valid{false};
		const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0},
		                _f{0};
	} SfcInteractionPar;

	/** Return a set of parameters which can be used to calculate the H2 formation rate on the
	    surface of grains of this type. */
	static SfcInteractionPar sfcInteractionPar(GasModule::GrainType t);

	/** Indicates wheter a photoelectric heating recipe has been implemented for the given type. */
	static bool heatingAvailable(GasModule::GrainType t);
};
#endif /* GASMODULE_GIT_SRC_GRAINTYPEPROPERTIES_H_ */
