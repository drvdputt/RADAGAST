#ifndef GASMODULE_GIT_SRC_SILICATEGRAIN_H_
#define GASMODULE_GIT_SRC_SILICATEGRAIN_H_

#include "GrainType.h"

class SilicateGrain : public GrainType
{
public:
	SilicateGrain();

	double photoElectricYield(double a, int z, double hnu) const override;

	double ionizationPotential(double a, int z) const override;

	double autoIonizationThreshold(double a) const override;

	double stickingCoefficient(double a, int z, int z_i) const override;
};

#endif /* GASMODULE_GIT_SRC_SILICATEGRAIN_H_ */
