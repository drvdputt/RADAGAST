#ifndef _SRC_GASINTERFACE_H_
#define _SRC_GASINTERFACE_H_

#define NO_SCATTER

#include "GasState.h"

#include <memory>

class GasInterfaceImpl;

/* The interface class that other codes should use. */
class GasInterface
{
public:
	GasInterface(const std::valarray<double>& frequencyv);
	GasInterface(const std::valarray<double>& frequencyv, bool improveGrid);
	~GasInterface();

	/* Exports the state of the gas as a compact, opaque object. Codes which make use of the gas
	 module can use these objects to store a number of gas states. They can repeatedly give
	 objects of this type to a single instance of the HydrogenCalculator to calculate any
	 optical properties on-the-fly. The exact contents of a GasState and the way the optical
	 properties are calculated (derived from densities vs caching them for example) are entirely
	 up to the implementations of the functions below and the definition in GasState.h. */
	void updateGasState(GasState& gs, double n, double Tinit,
	                    const std::valarray<double>& specificIntensityv) const;
	void initializeGasState(GasState& gs, double n, double T) const;

	/* The functions below hould provide a fast implementation to obtain the optical properties.
	 The implementation will depend on what is stored in a GasState object. A good balance
	 between the size of the GasState objects and the computation time needed for the optical
	 properties needs to be found. */
	// 1 erg / cm3 / s / Hz / sr = 0.1 J / m3 / s / Hz / sr
	double effectiveEmissivity_SI(const GasState& gs, size_t iFreq) const;
	double opacity_SI(const GasState& gs, size_t iFreq) const;
	double scatteringOpacity_SI(const GasState& gs, size_t iFreq) const;
	double absorptionOpacity_SI(const GasState& gs, size_t iFreq) const;

	std::valarray<double> frequencyv() const { return _frequencyv; }

	void testHeatingCurve(double n, const std::valarray<double>& specificIntensityv) const;

private:
	void zeroOpticalProperties(GasState& gs) const;

	std::valarray<double> _frequencyv;
	std::unique_ptr<GasInterfaceImpl> _pimpl;
};

#endif /* _SRC_GASINTERFACE_H_ */
