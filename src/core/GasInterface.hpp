#ifndef CORE_GASINTERFACE_HPP
#define CORE_GASINTERFACE_HPP

#include "GasState.hpp"
#include "GrainInterface.hpp"

#include <memory>
#include <string>

/** The real implementation is in GasInterfaceImpl(.h/.cpp). We only forward declare here, to
    hide as many dependencies as possible. Currently, a client code only needs GasState.h,
    GrainInterface.h, and this file in its include path. All public facing functionality (all
    classes and functions in the aforementioned files), falls under the namepace GasModule. */
class GasInterfaceImpl;

/** The full solution */
class GasSolution;

/** More details about the gas can be retrieved by passing an object of this class. Include its
    header to construct one, and then pass a pointer to the update function */
class GasDiagnostics;

namespace GasModule
{
/** The interface class that other codes should use. This class and GasState should be the only
    direct dependencies introduced in another code. Hence, the rest of the gas module code does
    not have to be recompiled whenever a change is made to a source file of the host code which
    makes use of this module. */
class GasInterface
{
public:
	// SETUP //

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

	GasInterface(GasInterface&&) = default;

	~GasInterface();

	std::valarray<double> iFrequencyv() const { return _iFrequencyv; }
	std::valarray<double> oFrequencyv() const { return _oFrequencyv; }
	std::valarray<double> eFrequencyv() const { return _eFrequencyv; }

	// USAGE 1: Calculate GasState objects //

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
	    needs to be tabulated for the frequencies contained in iFrequencyv. */
	void updateGasState(GasState&, double n,
	                    const std::valarray<double>& specificIntensityv,
	                    GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

	/** Does the same as the above, but without an input radiation field. Instead, a
	    blackbody of the given temperature is used to calculate GasState. It is recommended
	    to apply this function to all gas states before starting a simulation. */
	void initializeGasState(GasState&, double n, double T, GrainInterface&,
	                        GasDiagnostics* gd = nullptr) const;

	// USAGE 2: Translate GasState into optical properties //

	/** The functions below hould provide a fast implementation to obtain the optical
	    properties from the information stored in the GasState + constant information stored
	    in the GasInterface. */

	/** The emissivity in SI units, for a given frequency index. (converted using 1 erg cm-3
	    s-1 Hz-1 sr-1 = 0.1 J m-3 s-1 Hz-1 sr-1). [W m-3 Hz-1 sr-1] */
	double emissivity_SI(const GasState& gs, size_t iFreq) const;

	/** Returns the total opacity in SI units (converted from cm-1 = 100 * m-1). [m-1]*/
	double opacity_SI(const GasState& gs, size_t iFreq) const;

	// /** Returns the scattering opacity. Currently there is are no contributions, and this
	//     function simply returns zero for each prequency. [m-1] */
	// double scatteringOpacity_SI(const GasState& gs, size_t iFreq) const;

	// /** Return the absorption opacity. Currently, all opacity is absorption opacity and the
	//     albedo is therefore zero. [m-1] */
	// double absorptionOpacity_SI(const GasState& gs, size_t iFreq) const;

private:
	/** Modifies the contents of a GasState so that it is equivalent to no gas at all. */
	void zeroOpticalProperties(GasState& gs) const;

	std::valarray<double> _iFrequencyv;
	std::valarray<double> _oFrequencyv;
	std::valarray<double> _eFrequencyv;

	/* The implementation details, especially those that require the inclusion of other
	   files than 'GasState.h' in this header, are hidden behind this pointer. This way, the
	   dependency chain is broken off just before the top level of abstraction (i.e. this
	   header). This simplifies the include statements used for compiling a code which hosts
	   the gas module. */
	std::unique_ptr<GasInterfaceImpl> _pimpl;

public:
	/* This can be used to test individual components of the implementation. TODO: remove
	   this eventually.*/
	const GasInterfaceImpl* pimpl() const { return _pimpl.get(); }
};
} /* namespace GasModule */
#endif // CORE_GASINTERFACE_HPP
