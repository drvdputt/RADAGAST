#ifndef CORE_CARBONACEOUSGRAIN_H_
#define CORE_CARBONACEOUSGRAIN_H_

#include "GrainType.hpp"

/** This class provides properties for 'astonomical graphite' grains, and to some extent smaller
    carbonaceous grains such as some PAH's (at least, that's how it's done in the photoelectric
    heating recipe. */
class CarbonaceousGrain : public GrainType
{
public:
	CarbonaceousGrain();

	/** @name WD01 wrapper functions

	    These functions provide different properties used for the photoelectric heating recipe
	    of Weingartner \& Draine (2001). They are all wrappers around the functions defined in
	    \c WeingartnerDraine2001.h, calling the 'carbonaceous' versions of the latter by using
	    \c true for their \c carbonaceous argument. */
	/**@{*/
	double photoElectricYield(double a, int z, double hnu) const override;

	double ionizationPotential(double a, int z) const override;

	double autoIonizationThreshold(double a) const override;

	double stickingCoefficient(double a, int z, int z_i) const override;
	/**@}*/
};

#endif /* CORE_CARBONACEOUSGRAIN_H_ */
