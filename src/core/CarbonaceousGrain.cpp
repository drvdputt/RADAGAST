#include "CarbonaceousGrain.h"
#include "Constants.h"
#include "GasGrainInteraction.h"
#include "WeingartnerDraine2001.h"

CarbonaceousGrain::CarbonaceousGrain()
                : GrainType(GasModule::GrainTypeLabel::CAR, GasGrain::carSurface,
                            GasGrain::grainHeatingPerH2Formed_car, true,
                            WD01::workFunction(true))
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
