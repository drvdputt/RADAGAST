#include "GrainInterface.h"

#include "Array.h"
#include "Error.h"
#include "GrainTypeProperties.h"

namespace GasModule
{

GrainInterface::Population::Population(GrainType type, const Array& sizev, const Array& densityv,
                                       const Array& temperaturev, const std::vector<Array>& qAbsvv)
                : _type{type}, _sizev{sizev}, _densityv{densityv}, _temperaturev{temperaturev},
                  _qAbsvv{qAbsvv}, _h2FormationPars{GrainTypeProperties::sfcInteractionPar},
                  _heatingAvailable{GrainTypeProperties::heatingAvailable},
                  _workFunction{GrainTypeProperties::workFunction} _photoElectricYield{
                                  GrainTypeProperties::photoElectricYield},
                  _autoIonizationThreshold{GrainTypeProperties::autoIonizationThreshold},
                  _stickingCoefficient{GrainTypeProperties::stickingCoefficient}
{
	Error::equalCheck("sizev.size() and densityv.size()", sizev.size(), densityv.size());
	Error::equalCheck("sizev.size() and temperaturev.size()", sizev.size(),
	                  temperaturev.size());
	Error::equalCheck("sizev.size() and qAbsvv.size()", sizev.size(), qAbsvv.size());
}

GrainInterface::GrainInterface() = default;

GrainInterface::GrainInterface(const std::vector<Population>& populationv)
                : _populationv{populationv}
{
}

} /* namespace GasModule */
