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

GrainInterface::GrainInterface(std::unique_ptr<std::vector<GrainPopulation>> populationvToMove)
                : _populationv(std::move(populationvToMove))
{
}

size_t GrainInterface::numPopulations() const
{
	return _populationv ? _populationv->size() : 0;
}

const std::vector<GrainPopulation>* GrainInterface::populationv() const
{
	return _populationv.get();
}

void GrainInterface::test() const
{
	if (!_populationv)
		return;
	for (const auto& population : *_populationv)
		population.test();
}

} /* namespace GasModule */
