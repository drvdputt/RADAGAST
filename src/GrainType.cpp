#include "GrainType.h"
#include "CarbonaceousGrain.h"
#include "Constants.h"
#include "CustomGrain.h"
#include "Error.h"
#include "SilicateGrain.h"

GrainType::GrainType(const GasModule::SfcInteractionPar& sfcInteractionPar, bool heatingAvailable,
                     double workFunction)
                : _sfcInteractionPar(sfcInteractionPar), _heatingAvailable(heatingAvailable),
                  _workFunction(workFunction)
{
}

GrainType::~GrainType() = default;

std::unique_ptr<GrainType> GrainTypeFactory::makeBuiltin(GasModule::GrainTypeLabel t)
{
	if (t == GasModule::GrainTypeLabel::CAR)
		return std::make_unique<CarbonaceousGrain>();
	else if (t == GasModule::GrainTypeLabel::SIL)
		return std::make_unique<SilicateGrain>();
	else
	{
		Error::runtime("Only CAR and SIL have built-in grain properties.");
		return nullptr;
	}
}

std::unique_ptr<GrainType>
GrainTypeFactory::makeCustom(const GasModule::SfcInteractionPar& sfcInteractionPar,
                             bool heatingAvailable, double workFunction,
                             GasModule::IonizationPotentialf ionizationPotentialf,
                             GasModule::PhotoelectricYieldf photoelectricYieldf,
                             GasModule::AutoIonizationThresholdf autoIonizationThresholdf,
                             GasModule::StickingCoefficientf stickingCoefficientf)
{
	return std::make_unique<CustomGrain>(sfcInteractionPar, heatingAvailable, workFunction,
	                                     ionizationPotentialf, photoelectricYieldf,
	                                     autoIonizationThresholdf, stickingCoefficientf);
}
