#ifndef CORE_GASSTATE_HPP
#define CORE_GASSTATE_HPP

#include <string>
#include <valarray>

// Forward declarations for friend functions
namespace GasModule
{
class GasInterface;
class GasState;
} // namespace GasModule

class GasInterfaceImpl;

namespace Testing
{
void writeGasState(const std::string&, const GasModule::GasInterface&,
                   const GasModule::GasState&);
}

namespace GasModule
{
/** This class is part of the 'public' \c GasModule interface. It provides a way to store the
    output of the calculations done in the gas module in an opaque manner. A client code should
    typically have a list of \c GasState objects; one per computational volume element. The \c
    GasInterface is a friend class, as it is used to retrieve different properties from these
    gas states. By making these objects opaque, it is made sure that client codes do not depend
    on how the different properties are stored, and how the emissivity and opacity can be
    obtained. This way, we can change the way the result is stored and the spectra are
    calculated, which will help in finding a trade-off between storage space and computation
    time (e.g. storing whole spectra vs calculating them on the fly from a minimal set of
    properties). */
class GasState
{
	friend class GasInterface;
	friend void Testing::writeGasState(const std::string&, const GasModule::GasInterface&,
	                                   const GasModule::GasState&);

public:
	/** Creates an empty GasState object. The resulting object can't be used for anything,
	    except to allocate space, and handing it to the update function of the \c
	    GasInterface. */
	GasState() {}

	/** Constructor, subject to change */
	GasState(const std::valarray<double>& emissivityv,
	         const std::valarray<double>& opacityv, double T,
	         std::valarray<double>& densityv)
	                : _emissivityv(emissivityv), _opacityv(opacityv), _temperature(T),
	                  _densityv(densityv)
	{
	}
	// TODO use setters instead?

	/** Returns the gas temperature in K. */
	double temperature() const { return _temperature; }

private:
	/** Memory-heavy, but simple implementation: just store all the output */
	std::valarray<double> _emissivityv, _opacityv, _scatteringOpacityv;

	/** Some basic diagnostics. For more advanced diagnostics, import and use the
	    GasDiagnostics object. */
	double _temperature{0};
	std::valarray<double> _densityv{0};
};
} /* namespace GasModule */

#endif // CORE_GASSTATE_HPP
