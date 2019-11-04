#ifndef CORE_BUILTINGRAINTYPE_HPP
#define CORE_BUILTINGRAINTYPE_HPP

#include "GrainInterface.hpp"

/** Abstract class which described all the functions that should be provided, for a specific grain
    type. */
class GrainType
{
protected:
	/** Member values to be provided by the subclass constructor. */
	GrainType(GasModule::GrainTypeLabel l,
	          const GasModule::SfcInteractionPar& sfcInteractionPar, double heatPerH2,
	          bool heatingAvailable, double workFunction);

public:
	virtual ~GrainType();

	GasModule::GrainTypeLabel label() const { return _label; }

	GasModule::SfcInteractionPar sfcInteractionPar() const { return _sfcInteractionPar; }

	double heatPerH2() const { return _heatPerH2; };

	bool heatingAvailable() const { return _heatingAvailable; }

	double workFunction() const { return _workFunction; }

	virtual double ionizationPotential(double a, int Z) const = 0;

	virtual double photoelectricYield(double a, int z, double hnuDiff,
	                                  double Emin) const = 0;

	virtual double autoIonizationThreshold(double a) const = 0;

	virtual double stickingCoefficient(double a, int z, int z_i) const = 0;

private:
	GasModule::GrainTypeLabel _label;
	GasModule::SfcInteractionPar _sfcInteractionPar;
	double _heatPerH2;
	bool _heatingAvailable;
	double _workFunction;
};

/** Factory which can create one of the builtins, given a valid value for the enum. I used to
    have a method for custom grain types too, but it was getting too complicated to do such a
    thing in a general way. I recommend writing more ad-hoc code like the CarbonaceousGrain
    subclass. */
class GrainTypeFactory
{
public:
	/** Factory method for subclasses. */
	static std::unique_ptr<GrainType> makeBuiltin(GasModule::GrainTypeLabel t);
};

#endif // CORE_BUILTINGRAINTYPE_HPP
