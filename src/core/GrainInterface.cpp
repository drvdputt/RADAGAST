#include "GrainInterface.h"
#include "Array.h"
#include "Error.h"
#include "GrainType.h"

#include <cassert>

namespace GasModule
{

SfcInteractionPar::SfcInteractionPar(double EH2, double Es, double EHp, double EHc,
                                     double aSqrt, double nuH2, double nuHc, double F)
                : _valid{true}, _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt},
                  _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
{
}

GrainInterface::Population::Population(GrainTypeLabel type, const Array& sizev,
                                       const Array& densityv, const Array& temperaturev,
                                       const std::vector<Array>& qAbsvv)
                : _sizev{sizev}, _densityv{densityv},
                  _temperaturev{temperaturev}, _qAbsvv{qAbsvv},
                  _type(GrainTypeFactory::makeBuiltin(type))
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(),
	                  temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());
}

GrainInterface::Population::Population(const Array& sizev, const Array& densityv,
                                       const Array& temperaturev,
                                       const std::vector<Array>& qAbsvv,
                                       const SfcInteractionPar& sfcInteractionPar,
                                       bool heatingAvailable, double workFunction,
                                       IonizationPotentialf ionizationPotentialf,
                                       PhotoelectricYieldf photoElectricYieldf,
                                       AutoIonizationThresholdf autoIonizationThresholdf,
                                       StickingCoefficientf stickingCoefficientf)
                : _sizev{sizev}, _densityv{densityv},
                  _temperaturev{temperaturev}, _qAbsvv{qAbsvv},
                  _type(GrainTypeFactory::makeCustom(
                                  sfcInteractionPar, heatingAvailable, workFunction,
                                  ionizationPotentialf, photoElectricYieldf,
                                  autoIonizationThresholdf, stickingCoefficientf))
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(),
	                  temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());
}

GrainInterface::Population::Population(Population&& other) = default;

GrainInterface::Population::~Population() = default;

void GrainInterface::Population::test() const
{
	assert(_sizev.size() == _densityv.size());
	assert(_densityv.size() == _temperaturev.size());
	assert(_qAbsvv.size() == _sizev.size());
}

GrainInterface::GrainInterface() = default;

GrainInterface::~GrainInterface() = default;

GrainInterface::GrainInterface(GrainInterface&&) = default;

GrainInterface::GrainInterface(
                std::unique_ptr<std::vector<GrainInterface::Population>> populationvToMove)
                : _populationv(std::move(populationvToMove))
{
}

size_t GrainInterface::numPopulations() const
{
	return _populationv ? _populationv->size() : 0;
}

const GrainInterface::Population* GrainInterface::population(size_t i) const
{
	if (!_populationv)
		Error::runtime("No populations in this object!");
	if (i > _populationv->size())
		Error::runtime("Population index out of bounds.");
	return _populationv->data() + i;
}

const std::vector<GrainInterface::Population>* GrainInterface::populationv() const
{
	return _populationv.get();
}

void GrainInterface::test() const
{
	if (!_populationv)
		return;
	for (auto& population : *_populationv)
		population.test();
}

} /* namespace GasModule */
