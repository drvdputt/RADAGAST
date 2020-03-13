#include "GrainInterface.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "GrainPopulation.hpp"
#include "SpecialFunctions.hpp"
#include <cassert>

namespace GasModule
{
    GrainInterface::GrainInterface() = default;

    GrainInterface::~GrainInterface() = default;

    GrainInterface::GrainInterface(GrainInterface&&) = default;

    void GrainInterface::addPopulation(GasModule::GrainTypeLabel type, const Array& sizev, const Array& densityv,
                                       const Array& temperaturev, const Array& frequencyv,
                                       const std::vector<Array>& qAbsvv)
    {
        if (!_populationv) _populationv = std::make_unique<std::vector<GrainPopulation>>();
        _populationv->emplace_back(type, sizev, densityv, temperaturev, frequencyv, qAbsvv);
    }

    void GrainInterface::changePopulationDensityv(int p, const std::valarray<double>& densityv)
    {
        _populationv->at(p).setDensityv(densityv);
    }

    size_t GrainInterface::numPopulations() const { return _populationv ? _populationv->size() : 0; }

    const std::vector<GrainPopulation>* GrainInterface::populationv() const { return _populationv.get(); }

    void GrainInterface::test() const
    {
        if (!_populationv) return;
        for (const auto& population : *_populationv) population.test();
    }
}
