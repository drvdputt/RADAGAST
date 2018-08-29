#include "CustomGrain.h"

CustomGrain::CustomGrain(const GasModule::SfcInteractionPar& sfcInteractionPar,
                         bool heatingAvailable, double workFunction,
                         GasModule::IonizationPotentialf ionizationPotentialf,
                         GasModule::PhotoelectricYieldf photoElectricYieldf,
                         GasModule::AutoIonizationThresholdf autoIonizationThresholdf,
                         GasModule::StickingCoefficientf stickingCoefficientf)
                : GrainType(sfcInteractionPar, heatingAvailable, workFunction),
                  _ionizationPotentialf(ionizationPotentialf),
                  _photoElectricYieldf(photoElectricYieldf),
                  _autoIonizationThresholdf(autoIonizationThresholdf),
                  _stickingCoefficientf(stickingCoefficientf)
{
}

double CustomGrain::photoElectricYield(double a, int z, double hnu) const
{
	return _photoElectricYieldf(a, z, hnu);
}

double CustomGrain::ionizationPotential(double a, int z) const
{
	return _ionizationPotentialf(a, z);
}

double CustomGrain::autoIonizationThreshold(double a) const
{
	return _autoIonizationThresholdf(a);
}

double CustomGrain::stickingCoefficient(double a, int z, int z_i) const
{
	return _stickingCoefficientf(a, z, z_i);
}
