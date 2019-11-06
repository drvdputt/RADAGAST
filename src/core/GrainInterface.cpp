#include "GrainInterface.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "SpecialFunctions.hpp"

#include <cassert>

namespace GasModule
{

GrainInterface::GrainInterface()
                : _populationv(std::make_unique<std::vector<GrainPopulation>>())
{
}

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

const GrainPopulation* GrainInterface::population(size_t i) const
{
	if (!_populationv)
		Error::runtime("No populations in this object!");
	if (i > _populationv->size())
		Error::runtime("Population index out of bounds.");
	return _populationv->data() + i;
}

std::vector<GrainPopulation>* GrainInterface::populationv() const { return _populationv.get(); }

void GrainInterface::test() const
{
	if (!_populationv)
		return;
	for (auto& population : *_populationv)
		population.test();
}

} /* namespace GasModule */
