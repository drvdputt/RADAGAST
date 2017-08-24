#include "GrainInfo.h"
#include "Error.h"

namespace GasModule
{

GrainInfo::Population::Population(GrainType type, const Array& sizev, const Array& densityv,
                                  const Array& temperaturev, const std::vector<Array>& qAbsvv)
                : _type{type}, _sizev{sizev}, _densityv{densityv},
                  _temperaturev{temperaturev}, _qAbsvv{qAbsvv}
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(), temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());
}

GrainInfo::GrainInfo() = default;

GrainInfo::GrainInfo(const Array& carbonaceousGrainSizev, const Array& carbonaceousDensityv,
                     const std::vector<Array>& carbonaceousAbsQvv, const Array& carTv,
                     const Array& silicateGrainSizev, const Array& silicateDensityv,
                     const std::vector<Array>& silicateAbsQvv, const Array& silTv)
                : _populationv{{GrainType::CAR, carbonaceousGrainSizev, carbonaceousDensityv, carTv,
                                carbonaceousAbsQvv},
                               {GrainType::SIL, silicateGrainSizev, silicateDensityv, silTv,
                                silicateAbsQvv}}
{
}

} /* namespace GasModule */
