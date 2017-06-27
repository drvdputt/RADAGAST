#ifndef _SRC_GASINTERFACE_H_
#define _SRC_GASINTERFACE_H_

#include "GasState.h"

#include <memory>

class GasInterfaceImpl;

/* The interface class that other codes should use. This class and GasState should be the only
   direct dependencies introduced in another code. Hence, the rest of the gas module code does not
   have to be recompiled whenever a change is made to a source file of the host code which makes use
   of this module. */
class GasInterface
{
public:
	/* Creates an instance of the gas module. Since the module currently works using a frequency
	  grid, a sorted list (ascending) of frequencies to perform the calculations on should be
	  supplied. */
	GasInterface(const std::valarray<double>& frequencyv);
	/* A hacky implementation to make 'suggesting' frequency points possible. This constructor
	   uses the given frequency points, and adds some by itself before the frequency grid is
	   passed through to all the members. This probably needs to change. */
	GasInterface(const std::valarray<double>& frequencyv, bool improveGrid);
	~GasInterface();

	/* Exports the state of the gas as a compact, opaque object. Codes which make use of the gas
	   module can use these objects to store a number of gas states. They can then repeatedly
	   give objects of this type to a single instance of the HydrogenCalculator to calculate any
	   optical properties on-the-fly. The exact contents of a GasState and the way the optical
	   properties are calculated (derived from densities vs caching them for example) are
	   entirely up to the implementations of the functions below and the definition in
	   GasState.h. */
	void updateGasState(GasState& gs, double n, double Tinit,
	                    const std::valarray<double>& specificIntensityv) const;

	/* Does the same as the above, but without an input radiation field. Instead, a blackbody of
	   the given temperature is used to calculate GasState. It is recommended to apply this
	   function to all gas states before starting a simulation. */
	void initializeGasState(GasState& gs, double n, double T) const;

	/* The functions below hould provide a fast implementation to obtain the optical properties.
	   The implementation will depend on what is stored in a GasState object. A good balance
	   between the size of the GasState objects and the computation time needed for the optical
	   properties needs to be found. */

	/* Based on the information in the given GasState, this function returns the emissivity in
	   SI units, for a given frequency index. (converted from 1 erg / cm3 / s / Hz / sr = 0.1 J
	   / m3 / s / Hz / sr). */
	double emissivity_SI(const GasState& gs, size_t iFreq) const;

	/* Returns the total opacity in SI units (converted from 1 / cm = 100 * 1 / m). */
	double opacity_SI(const GasState& gs, size_t iFreq) const;

	/* Returns the scattering opacity. Currently there is are no contributions, and this
	   function simply returns zero for each prequency */
	double scatteringOpacity_SI(const GasState& gs, size_t iFreq) const;

	/* Return the absorption opacity. Currently, all opacity is absorption opacity and the
	   albedo is therefore zero. */
	double absorptionOpacity_SI(const GasState& gs, size_t iFreq) const;

	/* Returns the frequency grid used. This may come in handy later, for example for situations
	   where the code generates its own frequency points. */
	std::valarray<double> frequencyv() const { return _frequencyv; }

	/* Test function. Better move this somewhere else. */
	void testHeatingCurve(double n, const std::valarray<double>& specificIntensityv) const;

private:
	/* Modifies the contents of a GasState so that it is equivalent to no gas at all. */
	void zeroOpticalProperties(GasState& gs) const;

	/* The frequency grid used for the calculations, which is either given at construction, or
	   generated by the gas module code itself. */
	std::valarray<double> _frequencyv;

	/* The implementation details, especially those that require the inclusion of other files
	   than 'GasState.h' in this header, are hidden behind this pointer. This way, the
	   dependency chain is broken off just before the top level of abstraction (i.e. this
	   header). This simplifies the include statements used for compiling a code which hosts the
	   gas module. */
	std::unique_ptr<GasInterfaceImpl> _pimpl;
};

#endif /* _SRC_GASINTERFACE_H_ */
