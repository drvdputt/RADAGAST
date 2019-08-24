#include "SilicateGrain.h"
#include "Constants.h"
#include "GasGrainInteraction.h"
#include "WeingartnerDraine2001.h"

SilicateGrain::SilicateGrain()
                : GrainType(GasModule::GrainTypeLabel::SIL, GasGrain::silSurface,
                            GasGrain::grainHeatingPerH2Formed_sil, true,
                            WD01::workFunction(false))
{
}

double SilicateGrain::photoelectricYield(double a, int z, double hnuDiff, double Emin) const
{
	return WD01::yield(a, z, hnuDiff, Emin, false);
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
