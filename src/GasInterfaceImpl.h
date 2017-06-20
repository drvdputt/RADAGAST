#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.h"
#include "GasState.h"
#include "NLevel.h"

#include <memory>

class FreeBound;
class FreeFree;

/* This class contains the actual implementation of the Gas Interface, and currently is the backbone
   of the hole code. That might change once more than one element is included, as this class is
   currently a gathering of Hydrogen-related functions. It calculates the ionized fraction, level
   populations and net heating rate of hydrogen repeatedly, until the equilibrium temperature is
   found.

   Once new species are added, chemistry will come into play, and the hydrogen-specific parts of
   this class might be moved down in the call hierarchy. The top level should then consist of a
   calculation of the chemistry between the different species, followed by a loop which calculates
   the individual properties of the species. Then, the total net heating rate for the gas mixture
   should be calculated, taking into account chemical heating, as well as the heating/cooling
   contributions which the species provide independently, such as line cooling.

   For now, different parts of the code are grouped by the kind of process. When an instance of
   GasInterfaceImpl is constructed, a few member which provide implementations of the various
   processes are constructed too. Objects of the classes FreeFree, FreeBound and NLevel are created,
   which decribe respectively the Brehmsstrahlung, the recombination continuum and bound-bound
   processes (between energy levels of Hydrogen). During their construction, they each read the data
   they need to perform their functions.

   As such, creating an object of GasInterface constructs a GasInterfaceImpl instance, which in turn
   creates three other objects describing underlying processes. This initialization chain can be
   expanded, with for example the insertion of a chemistry layer between GasInterface and what is
   currently called GasInterfaceImpl (this name will probably change).

   Analogously to the initialization, most of the functions in this class simply delegate their
   calculation to one of the members Example: The emissivity function simply asks each member for
   the emissivity that is generated by the processes it models. The results are simply summed; it is
   the task of this class to keep track of the different contributions.

   Some members also have their own balance calculations (at the moment only the NLevel member).
   This way, the calculation of the level populations, is delegated to NLevel. When expanding the
   code, extra 'layers' in this approach could be inserted. Consider the hypothetical chemistry
   layer again. The chemistry layer will keep track of all the different species, and solve its
   chemistry network to obtain their abundances. It could then invoke the level balance calculations
   of each of the species individually. These species will then make use of their own NLevel member
   to calculate the statisitcal equilibrium. To obtain the emissivity, the chemistry layer w sum the
   emissivities of all the different species (e.g. H + H2), the species will sum the emissivity of
   their own processes (lines + continuum). */

class GasInterfaceImpl
{
public:
	typedef struct
	{
		double n, T, f;
		Array specificIntensityv;
		NLevel::Solution levelSolution;
	} Solution;

	/* Creates an object which can calculate the NLTE state, of a pure (ionized and atomic)
	   hydrogen gas and the resulting opacity and emission on the provided frequency grid. */
	GasInterfaceImpl(const Array& frequencyv);

	/* Creates a HydrogenCalculator and creates a frequency grid / adjusts the given
	   frequencygrid (by reference) if suggestedGrid = true. */
	GasInterfaceImpl(const Array& frequencyv, bool improveGrid);

	Array frequencyv() const { return _frequencyv; }

	~GasInterfaceImpl();

	/* Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature
	   guess, and a vector containing the radiation field in specific intensity per frequency
	   units (on the same frequency grid as the one provided at construction). */
	void solveBalance(GasState&, double n, double Tinit, const Array& specificIntensity) const;

	/* Used by the balance solver to calculate the ionization fraction and level populations for
	   a certain electron temperature, under influence of a blackbody isrf of that same
	   temperature. Can be used by the client to manually set the temperature and calculate some
	   properties which can be used as an initial guess. */
	void solveInitialGuess(GasState&, double n, double T) const;

private:
	/* Calculates all the densities for a fixed temperature. Is repeatedly called by this
	   class */
	Solution calculateDensities(double n, double T, const Array& specificIntensityv) const;

public:
	/* The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emissivityv(const Solution&) const;

	/* The total opacity at each frequency in 1 / cm */
	Array opacityv(const Solution&) const;

	/* A possible scattering opacity */
	Array scatteringOpacityv(const Solution&) const;

	/* The total bolometric emission, in erg / s / cm^3, obtained by integrating the
	   emissivity. */
	double cooling(const Solution&) const;

	/* The total bolometric absorption, in erg / s / cm^3. This is an integral of the opacity
	   times the radiation field. */
	double heating(const Solution&) const;

	/* The bolometric emission by the lines only. The emissivity of the photon re-emissions is
	   also included in this value. */
	double lineCooling(const Solution&) const;

	/* The bolometric absorption by the lines only. The absorption of re-emitted line photons is
	   also included here. */
	double lineHeating(const Solution&) const;

	/* The bolometric emission by the continuum only (= cooling by recombination + free-free
	   continuum) */
	double continuumCooling(const Solution&) const;

	/* The bolometric absorption by the continuum only (free-free heating) in erg / s / cm3 */
	double continuumHeating(const Solution&) const;

	inline double np_ne(const Solution& s) const
	{
		double np = s.n * s.f;
		return np * np;
	}

	inline double nAtomic(const Solution& s) const { return (1 - s.f) * s.n; }

	void testHeatingCurve(double n, const Array& specificIntensityv) const;

private:
	/* To be set in constructor */
	Array _frequencyv;

	/* Pointers to other parts of the implementation, to make late initialization possible */
	std::unique_ptr<NLevel> _boundBound;
	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound;
	std::unique_ptr<FreeFree> _freeFree;
};

#endif /* _HYDROGENCALCULATOR_H_ */
