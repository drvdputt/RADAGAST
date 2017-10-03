#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.h"
#include "GasState.h"
#include "GrainInterface.h"
#include "GrainPhotoelectricEffect.h"
#include "NLevel.h"

class ChemistrySolver;
class FreeBound;
class FreeFree;
class NLevel;
class H2Levels;

/** This class contains the actual implementation of the Gas Interface, and currently is the
    backbone of the whole code. That might change once more than one element is included, as this
    class is currently a gathering of Hydrogen-related functions. It calculates the ionized
    fraction, level populations and net heating rate of hydrogen repeatedly, until the equilibrium
    temperature is found.

    For now, different parts of the code are grouped by the kind of process. When an instance of
    GasInterfaceImpl is constructed, a few members which provide implementations of the various
    processes are constructed too. Objects of the classes FreeFree, FreeBound and NLevel are
    created, which decribe respectively the Brehmsstrahlung, the recombination continuum and
    bound-bound processes (between energy levels of Hydrogen). During their construction, they each
    read the data they need to perform their functions.

    The most recent addition is a chemistry network, which can calculate the abundances of the
    different species, give a set of rate coefficients for the different formation and destruction
    processes.  The invocations of the chemistry calculation, level calculation of H, and
    calculation of H2 are currently at the same level in the call hierarchy. The implementation
    alternates between chemistry and level population calculations, until a self-consistent result
    is reached. Functions which calculate the optical properties and thermal effects of H and H2 are
    explicitly invoked from this class.

    When more species which are optically or thermally active are added, the hydrogen-specific parts
    of this class might be moved down in the call hierarchy. The species should be intercompatible
    so they can be put in a list, and the functions which calculate their individual optical
    properties should be callable from a loop. The total net heating rate for the gas mixture should
    take into account chemical heating, as well as the heating/cooling contributions which the
    species provide independently. This clearly calls for polymorphism. The decision still has to be
    made on how we want to group the processes.  We could make a list of objects representing each
    species, and indicate which of the objects provides heating, cooling, line emission, continuum,
    emission, opacity... or we could make separate lists of of heating providers, line emission
    provides, opacity providers ... and do a multiple inheritance thing. Or a maybe even a lambda
    function approach, where we make lists of functions that need to be called to get the total
    heating, emission, opacity...

    Some members also have their own balance calculations (at the moment only the NLevel member).
    This way, the calculation of the level populations, is delegated to NLevel. It's pretty clear
    that we can put all the NLevel systems into a list, and call solveBalance. But different
    subclasses of NLevel provide some different information, such as HydrogenLevels for the
    two-photon continuum, or H2Levels for dissociation rates. */

class GasInterfaceImpl
{
public:
	typedef struct
	{
		Array specificIntensityv;
		double T;
		EVector speciesNv;
		NLevel::Solution HSolution;
		NLevel::Solution H2Solution;
	} Solution;

	/** Creates an object which can calculates the NLTE state, of a pure (ionized and atomic)
	    hydrogen gas and the resulting opacity and emission on the provided frequency grid. A
	    constructor for manual setup (an argument for each subprocess of which more than 1
	    usable subclass/configuration exists). The components can be set up outside of this
	    constructor, and ownership is then transferred using a unique pointer. */
	GasInterfaceImpl(std::unique_ptr<NLevel> hmodel, bool molecular, const Array& frequencyv);

	Array frequencyv() const { return _frequencyv; }

	~GasInterfaceImpl();

	/** Used by the balance solver to calculate the ionization fraction and level populations
	    for a certain electron temperature, under influence of a blackbody isrf of that same
	    temperature. Can be used by the client to manually set the temperature and calculate
	    some properties which can be used as an initial guess. */
	void solveInitialGuess(GasModule::GasState&, double n, double T,
	                       const GasModule::GrainInterface&) const;

	/** Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature
	    guess, and a vector containing the radiation field in specific intensity per frequency
	    units (on the same frequency grid as the one provided at construction). */
	void solveBalance(GasModule::GasState&, double n, double Tinit,
	                  const Array& specificIntensity, const GasModule::GrainInterface&) const;

	/** Calculates all the densities for a fixed temperature. Is repeatedly called by this class
	    until equilibrium is found. */
	Solution calculateDensities(double n, double T, const Array& specificIntensityv,
	                            const GasModule::GrainInterface&) const;

public:
	/** The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emissivityv(const Solution&) const;

	/** The total opacity at each frequency in 1 / cm */
	Array opacityv(const Solution&) const;

	/** A possible scattering opacity */
	Array scatteringOpacityv(const Solution&) const;

	/** The total bolometric emission, in erg / s / cm^3, obtained by integrating the
	    emissivity. */
	double cooling(const Solution&) const;

	/** Total heating sans grains */
	double heating(const Solution&) const;

	/** The total heating, including the grain photoelectric effect, in erg / s / cm^3. */
	double heating(const Solution&, const GasModule::GrainInterface&) const;

	/** The cooling by the lines only. */
	double lineCooling(const Solution&) const;

	/** The heating by the lines only. */
	double lineHeating(const Solution&) const;

	/** The cooling by recombination + free-free continuum */
	double continuumCooling(const Solution&) const;

	/** The heating by processes involving the continuum [erg / s / cm3]. */
	double continuumHeating(const Solution&) const;

	inline double np_ne(const Solution& s) const
	{
		return s.speciesNv(ine) * s.speciesNv(inp);
	}

	inline double f(const Solution& s) const
	{
		return s.speciesNv(inp) / (s.speciesNv(inH) + 2 * s.speciesNv(inH2));
	}

	inline double nAtomic(const Solution& s) const { return s.speciesNv[inH]; }

private:
	/* To be set in constructor */
	const Array& _frequencyv;

	// These are shorthand for ChemicalNetwork::speciesIndex.at["name"]
	int ine, inp, inH, inH2;
	std::unique_ptr<ChemistrySolver> _chemSolver;

	/* Pointers to other parts of the implementation, to make late initialization possible */
	std::unique_ptr<NLevel> _atomicLevels;
	std::unique_ptr<H2Levels> _molecular;
	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound;
	std::unique_ptr<FreeFree> _freeFree;
};

#endif /* _HYDROGENCALCULATOR_H_ */
