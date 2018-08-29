#ifndef GASMODULE_GIT_SRC_CUSTOMGRAIN_H_
#define GASMODULE_GIT_SRC_CUSTOMGRAIN_H_

#include "GrainType.h"

class CustomGrain : public GrainType
{
public:
	CustomGrain(const GasModule::SfcInteractionPar& sfcInteractionPar,
	            bool heatingAvailable, double workFunction,
	            GasModule::IonizationPotentialf ionizationpotentialf,
	            GasModule::PhotoelectricYieldf photoElectricYieldf,
	            GasModule::AutoIonizationThresholdf autoIonizationThresholdf,
	            GasModule::StickingCoefficientf stickingCoefficientf);

	double photoElectricYield(double a, int z, double hnu) const override;

	double ionizationPotential(double a, int z) const override;

	double autoIonizationThreshold(double a) const override;

	double stickingCoefficient(double a, int z, int z_i) const override;

private:
	// Photoelectric effect / grain charging
	GasModule::IonizationPotentialf _ionizationPotentialf;
	GasModule::PhotoelectricYieldf _photoElectricYieldf;
	GasModule::AutoIonizationThresholdf _autoIonizationThresholdf;
	GasModule::StickingCoefficientf _stickingCoefficientf;
};

#endif /* GASMODULE_GIT_SRC_CUSTOMGRAIN_H_ */
