#include "SpeciesModelManager.hpp"
#include "BigH2Model.hpp"
#include "HFromFiles.hpp"
#include "HHardCoded.hpp"
#include "Options.hpp"
#include "SimpleH2.hpp"
#include <sstream>

namespace GasModule
{
    SpeciesModelManager::SpeciesModelManager()
    {
        _hData = std::make_unique<HFromFiles>();

        if (Options::speciesmodelmanager_enableBigH2)
            _h2Data = std::make_unique<H2Data>(Options::speciesmodelmanager_bigH2_maxJ,
                                               Options::speciesmodelmanager_bigH2_maxV);
    }

    std::unique_ptr<HModel> SpeciesModelManager::makeHModel() const { return std::make_unique<HModel>(_hData.get()); }

    std::unique_ptr<H2Model> SpeciesModelManager::makeH2Model() const
    {
        if (Options::speciesmodelmanager_enableBigH2)
            return std::make_unique<BigH2Model>(_h2Data.get());
        else
            return std::make_unique<SimpleH2>(&_h2LTECool);
    }
}
