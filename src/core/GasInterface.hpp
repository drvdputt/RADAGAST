#ifndef CORE_GASINTERFACE_HPP
#define CORE_GASINTERFACE_HPP

#include "GasState.hpp"
#include "GrainInterface.hpp"

#include <memory>
#include <string>

class GasInterfaceImpl;
class GasSolution;
class GasDiagnostics;
class Spectrum;

namespace GasModule
{
/** The interface class that other codes should use. We use the PIMPL pattern here to minimize
    the amount of includes necessary on the client side, so that the internals of the gas module
    can be changed without having to recompile the client code. For the documentation, see
    GasInterfaceImpl. */
class GasInterface
{
public:
	/** Creates an instance of the gas module. Multiple frequency grids are used, of which 3
	    need to be specified by the user.

	    - iFrequencyv is the grid used to discretize the input radiation field (the specific
              intensity at a certain point in space, in [erg s-1 cm-1 Hz-1 units]).

	    - eFrequencyv will be the grid on which the emissivity is calculated.

	    - oFrequencyv is the grid that will be used to discretize the output opacity. This
	      is typically coarser because a radiative transfer algorithm usually needs the
	      opacity in each grid cell.

	    Some configuration options in the form of strings are also provided. These need to
	    be changed once the rest of the code is ready. */
	GasInterface(const std::valarray<double>& iFrequencyv,
	             const std::valarray<double>& oFrequencyv,
	             const std::valarray<double>& eFrequencyv,
	             const std::string& atomChoice = "",
	             const std::string& moleculeChoice = "");

	std::valarray<double> iFrequencyv() const;
	std::valarray<double> oFrequencyv() const;
	std::valarray<double> eFrequencyv() const;

	/** The main way to run the code for a cell. A minimal set of results is stored in the
	    given GasState object. The exact contents of the GasState are not known to the user.
	    Through this interface, the variable information contained in the gas state can be
	    combined with other constants and functions to retrieve the opacity and emissivity
	    of the gas.

	    The other arguments reflect the physical conditions in the cell for which a gas
	    state is being calculated. The density of hydrogen nuclei, n; the ambient radiation
	    field in [erg s-1 cm-1 Hz-1] units, as discretized on the @c iFrequencyv grid given
	    as an argument to the constructor, originally; and a GrainInterface object,
	    describing the properties of the grains within the cell for which the user wants to
	    solve the gas state. See the @c GrainInterface documentation for information on how
	    to construct one of these objects. The absorption efficiency of the grains currently
	    needs to be tabulated for the frequencies contained in iFrequencyv.

	    Note that the temperatures in GrainInterface can be modified, to take into account
	    the effect of gas-grain collisions. */
	void updateGasState(GasState&, double n,
	                    const std::valarray<double>& specificIntensityv,
	                    GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

	/** Does the same as the above, but without an input radiation field. Instead, a
	    blackbody of the given temperature is used to calculate GasState. It is recommended
	    to apply this function to all gas states before starting a simulation. */
	void initializeGasState(GasState&, double n, double T, GrainInterface&,
	                        GasDiagnostics* gd = nullptr) const;

	/** The emissivity in SI units, for a given frequency index. (converted using 1 erg cm-3
	    s-1 Hz-1 sr-1 = 0.1 J m-3 s-1 Hz-1 sr-1). [W m-3 Hz-1 sr-1] */
	double emissivity_SI(const GasState& gs, size_t iFreq) const;

	/** Returns the total opacity in SI units (converted from cm-1 = 100 * m-1). [m-1]*/
	double opacity_SI(const GasState& gs, size_t iFreq) const;

	GasSolution solveInitialGuess(double n, double T, GrainInterface&) const;

	GasSolution solveTemperature(double n, const Spectrum& specificIntensity,
	                             GasModule::GrainInterface&) const;

	GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity,
	                           GasModule::GrainInterface&,
	                           double h2FormationOverride = -1) const;

	void solveDensities(GasSolution&, double n, double T, const Spectrum& specificIntensity,
	                    GasModule::GrainInterface&, bool startFromCurrent = false,
	                    double h2FormationOverride = -1) const;

	GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& specificIntensity,
	                               GasModule::GrainInterface&) const;

private:
	/* The implementation details, especially those that require the inclusion of other
	   files than 'GasState.h' in this header, are hidden behind this pointer. This way, the
	   dependency chain is broken off just before the top level of abstraction (i.e. this
	   header). This simplifies the include statements used for compiling a code which hosts
	   the gas module. */
	std::unique_ptr<GasInterfaceImpl> _pimpl;
};
} // namespace GasModule
#endif // CORE_GASINTERFACE_HPP
