#include "GrainPopulation.hpp"
#include "Error.hpp"
#include "GrainH2Formation.hpp"
#include "GrainPhotoelectricData.hpp"

GrainPopulation::GrainPopulation(GrainTypeLabel type, const Array& sizev, const Array& densityv,
                                 const Array& temperaturev, const std::vector<Array>& qAbsvv)
                : _sizev{sizev}, _densityv{densityv},
                  _temperaturev{temperaturev}, _qAbsvv{qAbsvv}
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(),
	                  temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());

	if (type == GrainTypeLabel::CAR)
	{
		_h2formation = std::make_unique<GrainH2Formation>(
		                GrainH2FormationData::carSurface,
		                GrainH2FormationData::grainHeatingPerH2Formed_car);
		_photoelectricData = std::make_unique<GrainPhotoelectricData>(true);
	}
	else if (type == GrainTypeLabel::SIL)
	{
		_h2formation = std::make_unique<GrainH2Formation>(
		                GrainH2FormationData::silSurface,
		                GrainH2FormationData::grainHeatingPerH2Formed_sil);
		_photoelectricData = std::make_unique<GrainPhotoelectricData>(false);
	}
	// else, no data is available and these will be nullptr
}

void GrainPopulation::test() const
{
	assert(_sizev.size() == _densityv.size());
	assert(_densityv.size() == _temperaturev.size());
	assert(_qAbsvv.size() == _sizev.size());
	if (TemplatedUtils::contains(0., _sizev))
		Error::runtime("Grain of size 0 not allowed!");
}
