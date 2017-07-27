#include "GrainInfo.h"
#include "Error.h"

namespace GasModule
{
GrainInfo::GrainInfo() = default;

GrainInfo::GrainInfo(GrainType t, const Array& grainSizev, const std::vector<Array>& absQvv)
                : _hasCarbon(t == GrainType::CARBON), _hasSilica(t == GrainType::SILICA),
                  _carbonSizeDist{t == GrainType::CARBON ? GrainSizeDistribution(grainSizev, absQvv)
                                                 : GrainSizeDistribution{}},
                  _silicaSizeDist{t == GrainType::SILICA ? GrainSizeDistribution(grainSizev, absQvv)
                                                 : GrainSizeDistribution{}}
{
	Error::equalCheck("grainSizev.size() and absQvv.size()", grainSizev.size(), absQvv.size());
}

GrainInfo::GrainInfo(const Array& carbonGrainSizev, const std::vector<Array>& carbonAbsQvv,
                     const Array& silicaGrainSizev, const std::vector<Array>& silicaAbsQvv)
                : _hasCarbon{true}, _hasSilica{true}, _carbonSizeDist(carbonGrainSizev, carbonAbsQvv),
                  _silicaSizeDist(silicaGrainSizev, silicaAbsQvv)
{
	Error::equalCheck("carbonGrainSizev.size() and carbonAbsQvv.size()",
	                  carbonGrainSizev.size(), carbonAbsQvv.size());
	Error::equalCheck("silicaGrainSizev.size() and silicaAbsQvv.size()",
	                  silicaGrainSizev.size(), silicaAbsQvv.size());
}

} /* namespace GasModule */
