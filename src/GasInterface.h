#ifndef GASMODULE_GIT_SRC_GASINTERFACE_H_
#define GASMODULE_GIT_SRC_GASINTERFACE_H_

#include "GasState.h"
#include "GrainInterface.h"

#include <memory>
#include <string>

/** The real implementation is in GasInterfaceImpl(.h/.cpp). We only forward declare here, to
    hide as many dependencies as possible. Currently, a client code only needs GasState.h,
    GrainInterface.h, and this file in its include path. All public facing functionality (all
    classes and functions in the aforementioned files), falls under the namepace GasModule. */
class GasInterfaceImpl;

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

	    - frequencyv is used for all the rest, until we replace it by (maybe):

	    - eFrequencyv will be the grid on which the emissivity is calculated. Following a
              recent discussion on the future design of SKIRT, this can typically have a
              somewhat larger resolution, if you do your radiative transfer per cell for
              example.

	    - oFrequencyv is the grid that will be used to discretize the output opacity. This
	      is typically coarser because a radiative transfer algorithm usually needs the
	      opacity in each grid cell.

	    (other approaches are still considered)

	    Some configuration options in the form of strings are also provided. These are
	    subject to change, and it is currently best to look at the source code of this
	    function to see which options are available. */
	GasInterface(const std::valarray<double>& iFrequencyv,
		     const std::valarray<double>& frequencyv,
	             /* const std::valarray<double>& eFrequencyv, */
	             /* const std::valarray<double>& oFrequencyv, */
	             const std::string& atomChoice = "",
	             const std::string& moleculeChoice = "");

	GasInterface(GasInterface&&) = default;

	~GasInterface();

	/** Returns the frequency grid used. This may come in handy later, for example for
	    situations where the code generates its own frequency points. */
	std::valarray<double> frequencyv() const { return _frequencyv; }

	/** Return the frequency grid used for the specific intensity (the one the user gave at
	    construction. */
	std::valarray<double> iFrequencyv() const { return _iFrequencyv; }

	// USAGE 1: Calculate GasState objects //

	/** Exports the state of the gas as a compact, opaque object. Codes which make use of
	    the gas module can use these objects to store a number of gas states. They can then
	    repeatedly give objects of this type to a single instance of the HydrogenCalculator
	    to calculate any optical properties on-the-fly. The exact contents of a GasState and
	    the way the optical properties are calculated (derived from densities vs caching
	    them for example) are entirely up to the implementations of the functions below and
	    the definition in GasState.h.

	    The arguments to be provided here are a GasState object to which the results will be
	    written; the density of hydrogen nuclei, n; the ambient radiation field in [erg s-1
	    cm-1 Hz-1] units, as discretized on the @c iFrequencyv grid given as an argument to
	    the constructor, originally; and a GrainInterface object, describing the properties
	    of the grains within the cell for which the user wants to solve the gas state. See
	    the @c GrainInterface documentation for information on how to construct one of these
	    objects. The frequencies used to tabulate the absorption efficiency of the grains in
	    its constructor need to be the same as those given for the radiation field here. */
	void updateGasState(GasState&, double n, double Tinit,
	                    const std::valarray<double>& specificIntensityv,
	                    const GrainInterface&) const;

	/** Does the same as the above, but without an input radiation field. Instead, a
	    blackbody of the given temperature is used to calculate GasState. It is recommended
	    to apply this function to all gas states before starting a simulation. */
	void initializeGasState(GasState&, double n, double T, const GrainInterface&) const;

	// USAGE 2: Translate GasState into optical properties //

	/** The functions below hould provide a fast implementation to obtain the optical
	    properties. The implementation will depend on what is stored in a GasState object. A
	    good balance between the size of the GasState objects and the computation time
	    needed for the optical properties needs to be found. */

	/** Based on the information in the given GasState, this function returns the emissivity
	    in SI units, for a given frequency index. (converted from 1 erg / cm3 / s / Hz / sr
	    = 0.1 J / m3 / s / Hz / sr). [W m-3 Hz-1 sr-1] */
	double emissivity_SI(const GasState& gs, size_t iFreq) const;

	/** Returns the total opacity in SI units (converted from 1 / cm = 100 * 1 / m). [m-1]*/
	double opacity_SI(const GasState& gs, size_t iFreq) const;

	/** Returns the scattering opacity. Currently there is are no contributions, and this
	    function simply returns zero for each prequency. [m-1] */
	double scatteringOpacity_SI(const GasState& gs, size_t iFreq) const;

	/** Return the absorption opacity. Currently, all opacity is absorption opacity and the
	    albedo is therefore zero. [m-1] */
	double absorptionOpacity_SI(const GasState& gs, size_t iFreq) const;

private:
	/** Modifies the contents of a GasState so that it is equivalent to no gas at all. */
	void zeroOpticalProperties(GasState& gs) const;

	std::valarray<double> _iFrequencyv;
	std::valarray<double> _frequencyv;

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
#endif /* GASMODULE_GIT_SRC_GASINTERFACE_H_ */
