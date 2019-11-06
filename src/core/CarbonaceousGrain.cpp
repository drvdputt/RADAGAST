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
