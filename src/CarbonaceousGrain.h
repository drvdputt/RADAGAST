#ifndef GASMODULE_GIT_SRC_CARBONACEOUSGRAIN_H_
#define GASMODULE_GIT_SRC_CARBONACEOUSGRAIN_H_

#include "GrainType.h"

class CarbonaceousGrain : public GrainType
{
public:
	CarbonaceousGrain();
	
	double photoElectricYield(double a, int z, double hnu) const override;

	double ionizationPotential(double a, int z) const override;

	double autoIonizationThreshold(double a) const override;

	double stickingCoefficient(double a, int z, int z_i) const override;
};

#endif /* GASMODULE_GIT_SRC_CARBONACEOUSGRAIN_H_ */
