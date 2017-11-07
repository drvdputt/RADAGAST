#include "SilicateGrain.h"
#include "Constants.h"
#include "WeingartnerDraine2001.h"
// TODO: Implement everything, and maybe use file-local constants instead of members of the parent
// class. At least put the constants next to their names, and then hand them over.

SilicateGrain::SilicateGrain()
                : GrainType({320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10}, true,
                            8 / Constant::ERG_EV)
{
}

double SilicateGrain::photoElectricYield(double a, int z, double hnu) const
{
	return WD01::yield(a, z, hnu, false);
}

double SilicateGrain::ionizationPotential(double a, int z) const
{
	return WD01::ionizationPotential(a, z, false);
}

double SilicateGrain::autoIonizationThreshold(double a) const
{
	return WD01::autoIonizationThreshold(a, false);
}

double SilicateGrain::stickingCoefficient(double a, int z, int z_i) const
{
	return WD01::stickingCoefficient(a, z, z_i, false);
}
