#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include <vector>
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
		return _ionizedFraction;
	}

private:
// Private constructor, only to be used by friended class which acts as a factory and can fill in all the members
	GasState(const std::vector<double>& frequencyv, const std::vector<double>& emissivityv,
			const std::vector<double>& opacityv, const std::vector<double>& scatteringOpacityv,
			double T, double f) :
			_frequencyv(frequencyv), _emissivityv(emissivityv), _opacityv(opacityv), _scatteringOpacityv(
					scatteringOpacityv), _T(T), _ionizedFraction(f)
	{
	}

private:
// Memory-heavy, but simple implementation: just store all the output
	std::vector<double> _frequencyv, _emissivityv, _opacityv, _scatteringOpacityv;

// Some diagnostics which are publicly available
	double _T{0};
	double _ionizedFraction{0};
};

#endif /* _GASSTATE_H_ */
