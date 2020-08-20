#include "SpeciesModelManager.hpp"
#include "BigH2Model.hpp"
#include "HFromFiles.hpp"
#include "HHardCoded.hpp"
#include "Options.hpp"
#include "SimpleH2.hpp"
#include <sstream>

namespace RADAGAST
{
    SpeciesModelManager::SpeciesModelManager()
    {
        _hData = std::make_unique<HFromFiles>();

        if (Options::speciesmodelmanager_enableBigH2)
            _h2Data = std::make_unique<H2Data>(Options::h2data_X_maxJ, Options::h2data_X_maxV, Options::h2data_E_maxJ,
                                               Options::h2data_E_minV, Options::h2data_E_maxV);
    }

    std::unique_ptr<HModel> SpeciesModelManager::makeHModel(const Spectrum* meanIntensity) const
    {
        return std::make_unique<HModel>(_hData.get(), meanIntensity);
    }

    std::unique_ptr<H2Model> SpeciesModelManager::makeH2Model(const Spectrum* meanIntensity) const
    {
        if (Options::speciesmodelmanager_enableBigH2)
            return std::make_unique<BigH2Model>(_h2Data.get(), meanIntensity);
        else
            return std::make_unique<SimpleH2>(&_h2LTECool, meanIntensity);
    }
}
