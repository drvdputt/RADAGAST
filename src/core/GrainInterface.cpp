#include "GrainInterface.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "GrainType.hpp"
#include "SpecialFunctions.hpp"

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

GrainInterface::Population::Population(Population&&) = default;

GrainInterface::Population::~Population() = default;

void GrainInterface::Population::test() const
{
	assert(_sizev.size() == _densityv.size());
	assert(_densityv.size() == _temperaturev.size());
	assert(_qAbsvv.size() == _sizev.size());
	if (TemplatedUtils::contains(0., _sizev))
		Error::runtime("Grain of size 0 not allowed!");
}

void GrainInterface::Population::calculateTemperature(std::valarray<double> frequencyv,
                                                      std::valarray<double> specificIntensityv,
                                                      std::valarray<double> otherGrainHeat)
{
	for (size_t i = 0; i < _sizev.size(); i++)
	{
		double cross = Constant::PI * _sizev[i] * _sizev[i];
		double absorption = cross * TemplatedUtils::integrate<double, Array, Array>(
		                                            frequencyv,
		                                            _qAbsvv[i] * specificIntensityv);

		auto heating = [&](double T) -> int {
			Array blackbodyIntegrandv(frequencyv.size());
			for (size_t j = 0; j < frequencyv.size(); j++)
				blackbodyIntegrandv[j] =
				                _qAbsvv[i][j] *
				                SpecialFunctions::planck(frequencyv[j], T);

			double bbEmission = 0;
			if (T > 0.)
				bbEmission = cross *
				             TemplatedUtils::integrate<double, Array, Array>(
				                             frequencyv, blackbodyIntegrandv);
			if (bbEmission < absorption + otherGrainHeat[i])
				return 1;
			if (bbEmission > absorption + otherGrainHeat[i])
				return -1;
			else
				return 0;
		};
		_temperaturev[i] = TemplatedUtils::binaryIntervalSearch<double>(heating, 30.,
		                                                                1.e-3, 300, 1.);
		DEBUG("New temp for grain " << i << " " << _temperaturev[i] << " K\n");
	}
}

GrainInterface::GrainInterface() : _populationv(std::make_unique<std::vector<Population>>()) {}

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

std::vector<GrainInterface::Population>* GrainInterface::populationv() const
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
