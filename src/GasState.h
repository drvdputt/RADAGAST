#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include <string>
#include <valarray>

namespace GasModule
{
class GasInterface;
}

class GasInterfaceImpl;

namespace Testing
{
void runGasInterfaceImpl(const GasModule::GasInterface&, const std::string&, double, double, double,
                         double);
}

namespace GasModule
{
class GasState
{
	friend class GasInterface;
	friend class ::GasInterfaceImpl;
	friend void Testing::runGasInterfaceImpl(const GasModule::GasInterface&, const std::string&,
	                                         double, double, double, double);

public:
	GasState() {}

	/** Returns the gas temperature in K. */
	double temperature() const { return _temperature; }
	/** Returns the formation rate of H2, in cm-3 s-1 */
	double h2form() const { return _h2form; }
	/** Return the heating rate by the grains in erg cm-3 s-1 */
	double grainHeat() const { return _grainHeat; }

	/** Returns the density of a species included in the model. For
	    the moment, some "inside knowledge" about the indices is required (0:e-, 1:H+, 2:H,
	    3:H2), but of course a more human readable format is desired for the future. TODO: use
	    some kind of map instead of "secret" indices. */
	double density(int i) const { return _densityv[i]; }
	double density_SI(int i) const { return 1e6 * _densityv[i]; }

	double ionizedFraction() const
	{
		return _densityv[1] / (_densityv[1] + _densityv[2] /*+ 2 * _densityv[3]*/);
	}

private:
	/** Private constructor, only to be used by friended class which acts as a factory and can
	    fill in all the members. */
	GasState(const std::valarray<double>& previousISRFv,
	         const std::valarray<double>& emissivityv, const std::valarray<double>& opacityv,
	         const std::valarray<double>& scatteringOpacityv, double T,
	         const std::valarray<double>& densityv, double h2form, double grainHeat)
	                : _previousISRFv(previousISRFv), _emissivityv(emissivityv),
	                  _opacityv(opacityv), _scatteringOpacityv(scatteringOpacityv),
		_temperature(T), _densityv(densityv), _h2form(h2form), _grainHeat(grainHeat)
	{
	}

private:
	/** Memory-heavy, but simple implementation: just store all the output */
	std::valarray<double> _previousISRFv, _emissivityv, _opacityv, _scatteringOpacityv;

	/** Some diagnostics which are publicly available */
	double _temperature{0};
	std::valarray<double> _densityv{0};
	double _h2form{0}, _grainHeat{0};
};
} /* namespace GasModule */

#endif /* _GASSTATE_H_ */
