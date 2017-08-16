#include "GrainInfo.h"
#include "Error.h"

namespace GasModule
{
GrainInfo::GrainInfo() = default;

GrainInfo::GrainInfo(GrainType t, const Array& grainSizev, const Array& grainDensityv,
                     const std::vector<Array>& absQvv, const Array& Tv)
                : _hasCar{t == GrainType::CAR}, _hasSil{t == GrainType::SIL},
                  _carSizeDist{t == GrainType::CAR
                                               ? GrainSizeDistribution(
                                                                 grainSizev, grainDensityv, absQvv,
                                                                 Tv,
                                                                 sfcInteractionPar(GrainType::CAR))
                                               : GrainSizeDistribution{}},
                  _silSizeDist{t == GrainType::SIL
                                               ? GrainSizeDistribution(
                                                                 grainSizev, grainDensityv, absQvv,
                                                                 Tv,
                                                                 sfcInteractionPar(GrainType::SIL))
                                               : GrainSizeDistribution{}}
{
	Error::equalCheck("grainSizev.size() and absQvv.size()", grainSizev.size(), absQvv.size());
	Error::equalCheck("grainSizev.size() and grainDensityv.size()", grainSizev.size(),
	                  grainDensityv.size());
}

GrainInfo::GrainInfo(const Array& carbonaceousGrainSizev, const Array& carbonaceousDensityv,
                     const std::vector<Array>& carbonaceousAbsQvv, const Array& carTv,
                     const Array& silicateGrainSizev, const Array& silicateDensityv,
                     const std::vector<Array>& silicateAbsQvv, const Array& silTv)
                : _hasCar{true}, _hasSil{true},
                  _carSizeDist(carbonaceousGrainSizev, carbonaceousDensityv, carbonaceousAbsQvv,
                               carTv, sfcInteractionPar(GrainType::CAR)),
                  _silSizeDist(silicateGrainSizev, silicateDensityv, silicateAbsQvv, silTv,
                               sfcInteractionPar(GrainType::SIL))
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
