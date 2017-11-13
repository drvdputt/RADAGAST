#include "GrainInterface.h"
#include "Array.h"
#include "Error.h"
#include "GrainType.h"

namespace GasModule
{

SfcInteractionPar::SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt,
                                     double nuH2, double nuHc, double F)
                : _valid{true}, _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt},
                  _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
{
}

GrainInterface::Population::Population(GrainTypeLabel type, const Array& sizev,
                                       const Array& densityv, const Array& temperaturev,
                                       const std::vector<Array>& qAbsvv)
                : _sizev{sizev}, _densityv{densityv}, _temperaturev{temperaturev}, _qAbsvv{qAbsvv},
                  _type{GrainTypeFactory::makeBuiltin(type)}
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(),
	                  temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());
}

GrainInterface::Population::Population(Population&& other) = default;

GrainInterface::Population::~Population() = default;

GrainInterface::GrainInterface() = default;

GrainInterface::~GrainInterface() = default;

GrainInterface::GrainInterface(const std::vector<GrainInterface::Population>* populationv)
                : _populationv{populationv}
{
}

int GrainInterface::numPopulations() const { return _populationv ? _populationv->size() : 0; }

const GrainInterface::Population* GrainInterface::population(int i) const
{
	if (!_populationv)
		Error::runtime("No populations in this object!");
	return &(*_populationv)[i];
}

const std::vector<GrainInterface::Population>* GrainInterface::populationv() const
{
	return _populationv;
}

} /* namespace GasModule */
