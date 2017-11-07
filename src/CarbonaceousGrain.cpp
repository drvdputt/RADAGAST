#include "CarbonaceousGrain.h"
#include "Constants.h"
#include "WeinGartnerDraine2001.h"

CarbonaceousGrain::CarbonaceousGrain()
                : GrainType({520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10}, true,
                            4.4 / Constant::ERG_EV)
{
}

double CarbonaceousGrain::photoElectricYield(double a, int z, double hnu) const
{
	return WD01::yield(a, z, hnu, true);
}
double CarbonaceousGrain::ionizationPotential(double a, int z) const
{
	return WD01::ionizationPotential(a, z, true);
}
double CarbonaceousGrain::autoIonizationThreshold(double a) const
{
	return WD01::autoIonizationThreshold(a, true);
}
double CarbonaceousGrain::stickingCoefficient(double a, int z, int z_i) const
{
	return WD01::stickingCoefficient(a, z, z_i, true);
}
