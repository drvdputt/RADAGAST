#ifndef _SRC_GASINTERFACE_H_
#define _SRC_GASINTERFACE_H_

#include "GasState.h"

#include <memory>

class HydrogenCalculator;

/* The interface class that other codes should use. */
class GasInterface
{
public:
	GasInterface(const std::valarray<double>& frequencyv);

	~GasInterface();

	/* Exports the state of the gas as a compact, opaque object. Codes which make use of the gas
	 module can use these objects to store a number of gas states. They can repeatedly give
	 objects of this type to a single instance of the HydrogenCalculator to calculate any
	 optical properties on-the-fly. The exact contents of a GasState and the way the optical
	 properties are calculated (derived from densities vs caching them for example) are entirely
	 up to the implementations of the functions below and the definition in GasState.h. */
	void updateGasState(GasState& gs, double n, double Tinit,
	                    const std::valarray<double>& specificIntensityv);
	void initializeGasState(GasState& gs, double n, double T);

	/* The functions below hould provide a fast implementation to obtain the optical properties.
	 The implementation will depend on what is stored in a GasState object. A good balance
	 between the size of the GasState objects and the computation time needed for the optical
	 properties needs to be found. */
	// 1 erg / cm3 = 0.1 J / m3
	double effectiveEmissivity_SI(const GasState& gs, size_t iFreq) const
	{
		double r = 0.1 * (gs._emissivityv[iFreq] /*-
		                  gs._scatteringOpacityv[iFreq] * gs._previousISRFv[iFreq]*/);
		return r > 0 ? r : 0;
	}

	// 1 / cm = 100 / m
	double opacity_SI(const GasState& gs, size_t iFreq) const
	{
		return 100 * gs._opacityv[iFreq];
	}
	double scatteringOpacity_SI(const GasState& gs, size_t iFreq) const
	{
		return 0;//100 * gs._scatteringOpacityv[iFreq];
	}
	double absorptionOpacity_SI(const GasState& gs, size_t iFreq) const
	{
		return 100 * (gs._opacityv[iFreq] /*- gs._scatteringOpacityv[iFreq]*/);
	}

private:
	void zeroOpticalProperties(GasState& gs) const;

	std::valarray<double> _frequencyv;
	std::unique_ptr<HydrogenCalculator> _hc{nullptr};
};

#endif /* _SRC_GASINTERFACE_H_ */
