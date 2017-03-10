#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include <vector>

class GasState
{
	friend class HydrogenCalculator;
private:
	// Private constructor, only to be used by friended class which acts as a factory and can fill in all the members
	GasState(const std::vector<double>& frequencyv, const std::vector<double>& emissivityv,
			const std::vector<double>& opacityv, const std::vector<double>& scatteringOpacityv)
	: _frequencyv(frequencyv), _emissivityv(emissivityv), _opacityv(opacityv), _scatteringOpacityv(scatteringOpacityv)
	{}

	// Memory-heavy, but simple implementation: just store all the output
	std::vector<double> _frequencyv, _emissivityv, _opacityv, _scatteringOpacityv;
};

#endif /* _GASSTATE_H_ */
