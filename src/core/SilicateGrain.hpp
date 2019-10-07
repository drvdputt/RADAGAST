#ifndef GASMODULE_GIT_SRC_SILICATEGRAIN_H_
#define GASMODULE_GIT_SRC_SILICATEGRAIN_H_

#include "GrainType.hpp"

/** This class provides properties for 'astronomical silicate' grains. */
class SilicateGrain : public GrainType
{
public:
	SilicateGrain();

	/** @name WD01 wrapper functions

	    These functions provide different properties used for the photoelectric heating recipe
	    of Weingartner \& Draine (2001). They are all wrappers around the functions defined in
	    \c WeingartnerDraine2001.h, calling the 'silicate' versions of the latter by using \c
	    false for their \c carbonaceous argument. */
	/**@{*/
	double photoElectricYield(double a, int z, double hnu) const override;

	double ionizationPotential(double a, int z) const override;

	double autoIonizationThreshold(double a) const override;

	double stickingCoefficient(double a, int z, int z_i) const override;
	/**@}*/
};

#endif /* GASMODULE_GIT_SRC_SILICATEGRAIN_H_ */
