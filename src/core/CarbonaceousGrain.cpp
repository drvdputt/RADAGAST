#include "CarbonaceousGrain.hpp"
#include "Constants.hpp"
#include "GasGrainInteraction.hpp"
#include "WeingartnerDraine2001.hpp"

CarbonaceousGrain::CarbonaceousGrain()
                : GrainType(GasModule::GrainTypeLabel::CAR, GasGrain::carSurface,
                            GasGrain::grainHeatingPerH2Formed_car, true,
                            WD01::workFunction(true))
{
}

double CarbonaceousGrain::photoelectricYield(double a, int z, double hnuDiff, double Emin) const
{
	return WD01::yield(a, z, hnuDiff, Emin, true);
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
