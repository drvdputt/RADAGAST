#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include "Array.h"
#include <iostream>

class GasState
{
	friend class HydrogenCalculator;

public:
	GasState()
	{
	}

	double temperature()
	{
		return _T;
	}
	double ionizedFraction()
	{
		return _f;
	}

private:
	/* Private constructor, only to be used by friended class which acts as a factory and can fill in all
	 the members. */
	GasState(const Array& frequencyv, const Array& previousISRFv, const Array& emissivityv,
			const Array& opacityv, const Array& scatteringOpacityv, double T, double f) :
			_frequencyv(frequencyv), _previousISRFv(previousISRFv), _emissivityv(emissivityv), _opacityv(
					opacityv), _scatteringOpacityv(scatteringOpacityv), _T(T), _f(f)
	{
	}

private:
	/* Memory-heavy, but simple implementation: just store all the output */
	Array _frequencyv, _previousISRFv, _emissivityv, _opacityv, _scatteringOpacityv;

	/* Some diagnostics which are publicly available */
	double _T
	{ 0 };
	double _f
	{ 0 };
};

#endif /* _GASSTATE_H_ */
