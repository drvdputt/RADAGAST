#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include <string>
#include <valarray>

class GasInterface;

namespace Testing
{
void runGasInterfaceImpl(const GasInterface&, const std::string&, double, double, double, double);
}

class GasState
{
	friend class GasInterface;
	friend class GasInterfaceImpl;
	friend void Testing::runGasInterfaceImpl(const GasInterface&, const std::string&, double,
	                                         double, double, double);

public:
	GasState() {}

	double temperature() { return _temperature; }
	double ionizedFraction() { return _ionizedFraction; }

private:
	/** Private constructor, only to be used by friended class which acts as a factory and can
	    fill in all the members. */
	GasState(const std::valarray<double>& previousISRFv,
	         const std::valarray<double>& emissivityv, const std::valarray<double>& opacityv,
	         const std::valarray<double>& scatteringOpacityv, double T, double f)
	                : _previousISRFv(previousISRFv), _emissivityv(emissivityv),
	                  _opacityv(opacityv), _scatteringOpacityv(scatteringOpacityv),
	                  _temperature(T), _ionizedFraction(f)
	{
	}

private:
	/** Memory-heavy, but simple implementation: just store all the output */
	std::valarray<double> _previousISRFv, _emissivityv, _opacityv, _scatteringOpacityv;

	/** Some diagnostics which are publicly available */
	double _temperature{0};
	double _ionizedFraction{0};
};

#endif /* _GASSTATE_H_ */
