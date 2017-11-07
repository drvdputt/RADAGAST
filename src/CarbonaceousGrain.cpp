#include "CarbonaceousGrain.h"
#include "Constants.h"

CarbonaceousGrain::CarbonaceousGrain()
                : GrainType({520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10}, true,
                            4.4 / Constant::ERG_EV)
{
}

double CarbonaceousGrain::photoElectricYield(double a, int z, double hnu) const { return 0.; }

double CarbonaceousGrain::ionizationPotential(double a, int z) const { return 0.; }

double CarbonaceousGrain::autoIonizationThreshold(double a) const { return 0; }

double CarbonaceousGrain::stickingCoefficient(double a, int z, int z_i) const { return 0; }
