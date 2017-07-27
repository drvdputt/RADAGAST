#include "GrainInfo.h"

namespace GasModule
{
GrainInfo::GrainInfo() = default;

GrainInfo::GrainInfo(GrainType t, const std::vector<double>& grainSizev,
                     const std::vector<std::vector<double>>& absQvv)
                : _contents{SingleTypeAvailable(t)}, _carbon{t == GrainType::CARBON ? GrainData{grainSizev, absQvv}
                                                               : GrainData{}},
                  _silica{t == GrainType::SILICA ? GrainData{grainSizev, absQvv} : GrainData{}}
{
}

GrainInfo::GrainInfo(const std::vector<double>& carbonGrainSizev,
                     const std::vector<std::vector<double>>& carbonAbsQvv,
                     const std::vector<double>& silicaGrainSizev,
                     const std::vector<std::vector<double>>& silicaAbsQvv)
                : _contents{AvailableContents::BOTH}, _carbon(carbonGrainSizev, carbonAbsQvv),
                  _silica(silicaGrainSizev, silicaAbsQvv)
{
}

} /* namespace GasModule */
