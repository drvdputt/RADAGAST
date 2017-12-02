#ifndef _GASSTATE_H_
#define _GASSTATE_H_

#include <string>
#include <valarray>

// Forward declarations for friend functions
namespace GasModule
{
class GasInterface;
class GasState;
}

class GasInterfaceImpl;

namespace Testing
{
void writeGasState(const std::string&, const GasModule::GasInterface&,
                   const GasModule::GasState&);
}

namespace GasModule
{

/** This class is part of the 'public' \c GasModule interface. It provides a way to store the output
    of the calculations done in the gas module in an opaque manner. A client code should typically
    have a list of \c GasState objects; one per computational volume element. The \c GasInterface is
    a friend class, as it is used to retrieve different properties from these gas states. By making
    these objects opaque, it is made sure that client codes do not depend on how the different
    properties are stored, and how the emissivity and opacity can be obtained. This way, we can
    change the way the result is stored and the spectra are calculated, which will help in finding a
    trade-off between storage space and computation time (e.g. storing whole spectra vs calculating
    them on the fly from a minimal set of properties). */
class GasState
{
	friend class GasInterface;
	friend class ::GasInterfaceImpl;
	friend void Testing::writeGasState(const std::string&, const GasModule::GasInterface&,
	                                   const GasModule::GasState&);

public:
	/** Creates an empty GasState object. The resulting object can't be used for anything,
	    except to allocate space, and handing it to the update function of the \c GasInterface. */
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

	/** Return the ionized fraction (protons / (proton + neutral)). */
	double ionizedFraction() const
	{
		return _densityv[1] / (_densityv[1] + _densityv[2] /*+ 2 * _densityv[3]*/);
	}

private:
	/** Private constructor, only to be used by friended class which acts as a factory and can
	    fill in all the members. */
	GasState(const std::valarray<double>& previousISRFv,
	         const std::valarray<double>& emissivityv,
	         const std::valarray<double>& opacityv,
	         const std::valarray<double>& scatteringOpacityv, double T,
	         const std::valarray<double>& densityv, double h2form, double grainHeat)
	                : _previousISRFv(previousISRFv), _emissivityv(emissivityv),
	                  _opacityv(opacityv), _scatteringOpacityv(scatteringOpacityv),
	                  _temperature(T), _densityv(densityv), _h2form(h2form),
	                  _grainHeat(grainHeat)
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
