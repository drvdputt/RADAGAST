#ifndef GASMODULE_GIT_SRC_BUILTINGRAINTYPE_H_
#define GASMODULE_GIT_SRC_BUILTINGRAINTYPE_H_

#include "GrainInterface.h"

/** Abstract class which described all the functions that should be provided, for a specific grain
    type. */
class GrainType
{
protected:
	/** Member values to be provided by the subclass constructor. */
	GrainType(const GasModule::SfcInteractionPar& sfcInteractionPar, bool heatingAvailable,
	          double workFunction);

public:
	~GrainType();

	GasModule::SfcInteractionPar sfcInteractionPar() const { return _sfcInteractionPar; }

	bool heatingAvailable() const { return _heatingAvailable; }

	double workFunction() const { return _workFunction; }

	virtual double ionizationPotential(double a, int Z) const = 0;

	virtual double photoElectricYield(double a, int z, double hnu) const = 0;

	virtual double autoIonizationThreshold(double a) const = 0;

	virtual double stickingCoefficient(double a, int z, int z_i) const = 0;

private:
	GasModule::SfcInteractionPar _sfcInteractionPar;
	bool _heatingAvailable;
	double _workFunction;
};

/** Factory which can create either a custom grain type, or on of the builtins, given a valid value
    for the enum (this is enforced at compile-time though, thanks to the use of enum class. */
class GrainTypeFactory
{
public:
	/** Factory method for subclasses. */
	static std::unique_ptr<GrainType> makeBuiltin(GasModule::GrainTypeLabel t);

	static std::unique_ptr<GrainType>
	makeCustom(const GasModule::SfcInteractionPar& sfcInteractionPar, bool heatingAvailable,
	           double workFunction, GasModule::IonizationPotentialf ionizationPotentialf,
	           GasModule::PhotoelectricYieldf photoelectricYieldf,
	           GasModule::AutoIonizationThresholdf autoIonizationThresholdf,
	           GasModule::StickingCoefficientf stickingCoefficientf);
};

#endif /* GASMODULE_GIT_SRC_BUILTINGRAINTYPE_H_ */
