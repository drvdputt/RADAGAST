#include "GrainInfo.h"
#include "Error.h"

namespace GasModule
{
GrainInfo::GrainInfo() = default;

GrainInfo::GrainInfo(GrainType t, const Array& grainSizev, const Array& grainDensityv,
                     const std::vector<Array>& absQvv)
                : _hasCar{t == GrainType::CARBONACEOUS}, _hasSil{t == GrainType::SILICATE},
                  _carSizeDist{t == GrainType::CARBONACEOUS
                                               ? GrainSizeDistribution(grainSizev, grainDensityv,
                                                                       absQvv)
                                               : GrainSizeDistribution{}},
                  _silSizeDist{t == GrainType::SILICATE
                                               ? GrainSizeDistribution(grainSizev, grainDensityv,
                                                                       absQvv)
                                               : GrainSizeDistribution{}}
{
	Error::equalCheck("grainSizev.size() and absQvv.size()", grainSizev.size(), absQvv.size());
	Error::equalCheck("grainSizev.size() and grainDensityv.size()", grainSizev.size(),
	                  grainDensityv.size());
}

GrainInfo::GrainInfo(const Array& carbonaceousGrainSizev, const Array& carbonaceousDensityv,
                     const std::vector<Array>& carbonaceousAbsQvv, const Array& silicateGrainSizev,
                     const Array& silicateDensityv, const std::vector<Array>& silicateAbsQvv)
                : _hasCar{true}, _hasSil{true},
                  _carSizeDist(carbonaceousGrainSizev, carbonaceousDensityv, carbonaceousAbsQvv),
                  _silSizeDist(silicateGrainSizev, silicateDensityv, silicateAbsQvv)
{
	Error::equalCheck("carbonaceousGrainSizev.size() and carbonaceousAbsQvv.size()",
	                  carbonaceousGrainSizev.size(), carbonaceousAbsQvv.size());
	Error::equalCheck("carbonaceousGrainSizev.size() and carbonaceousDensityv.size()",
	                  carbonaceousGrainSizev.size(), carbonaceousDensityv.size());
	Error::equalCheck("silicateGrainSizev.size() and silicateAbsQvv.size()",
	                  silicateGrainSizev.size(), silicateAbsQvv.size());
	Error::equalCheck("silicateGrainSizev.size() and silicateDensityv.size()",
	                  silicateGrainSizev.size(), silicateDensityv.size());
}

} /* namespace GasModule */
