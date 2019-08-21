#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.h"
#include "GasState.h"
#include "GrainInterface.h"
#include "GrainPhotoelectricEffect.h"
#include "NLevel.h"
#include "Spectrum.h"

class ChemistrySolver;
class FreeBound;
class FreeFree;
class GasDiagnostics;
class HydrogenLevels;
class H2Levels;

/** This class contains the actual implementation of the Gas Interface, and currently is the
    backbone of the whole code. It calculates the abundances, level populations and net heating
    rate of hydrogen repeatedly, until the equilibrium temperature is found.

    For now, different parts of the code are grouped by the kind of process. When an instance of
    GasInterfaceImpl is constructed, a few members which provide implementations of the various
    processes are constructed too. Objects of the classes FreeFree, FreeBound and NLevel are
    created, which decribe respectively the Brehmsstrahlung, the recombination continuum and
    bound-bound processes (between energy levels of Hydrogen). During their construction, they
    each read the data they need to perform their functions.

    The most recent addition is a chemistry network, which can calculate the abundances of the
    different species, given a set of rate coefficients for the different formation and
    destruction processes. The invocations of the chemistry calculation, level population
    calculations for H2 and H, and deriving new rate coefficients are currently at the same
    level in the call hierarchy. The implementation alternates between chemistry and level
    population calculations, until a self-consistent result is reached. Functions which
    calculate the optical properties and thermal effects of H and H2 are explicitly invoked from
    this class (as opposed to, for example, putting all the "opacity-havers" in a list).

    The total net heating rate for the gas mixture should take into account chemical heating, as
    well as the heating/cooling contributions which the species provide independently.

    The decision still has to be made on how we want to group the processes. We could make a
    list of objects representing each species, and indicate which of the objects provides
    heating, cooling, line emission, continuum, emission, opacity... or we could make separate
    lists of of heating providers, line emission providers, opacity providers ... (this smells of 
    multiple inheritance).

    It's pretty clear that we can put all the NLevel systems into a list, and call solveBalance.
    So solveBalance should have the same signature for all subclasses of NLevel. But different
    subclasses of NLevel provide some different information, such as HydrogenLevels for the
    two-photon continuum, or H2Levels for dissociation rates. This explains the need for
    subclassing. */
class GasInterfaceImpl
{
public:
	/** The main routine returns a \c Solution object, mainly to be used locally within this
	    class. For a solution more suitable for passing around the gas properties in general
	    between functions, check out GasStruct (TODO: merge these two concepts, maybe,
	    preferably). This struct contains a bunch of information about the final state,
	    including similar \c Solution objects for the levels. By passing around this object,
	    the operation of the module can be kept thread-safe. */
	typedef struct Solution
	{
		Spectrum specificIntensity;
		double T;
		EVector speciesNv;
		NLevel::Solution HSolution;
		NLevel::Solution H2Solution;
	} Solution;

	/** Creates an object which can calculates the NLTE state, of a pure (ionized and
	    atomic) hydrogen gas and the resulting opacity and emission on the provided
	    frequency grid. A constructor for manual setup (an argument for each subprocess of
	    which more than 1 usable subclass/configuration exists). The components can be set
	    up outside of this constructor, and ownership is then transferred using a unique
	    pointer. */
	GasInterfaceImpl(std::unique_ptr<HydrogenLevels> atomModel,
	                 std::unique_ptr<H2Levels> molecularModel);

	~GasInterfaceImpl();

	/** Used by the balance solver to calculate the ionization fraction and level
	    populations for a certain electron temperature, under influence of a blackbody isrf
	    of that same temperature. Can be used by the client to manually set the temperature
	    and calculate some properties which can be used as an initial guess. */
	Solution solveInitialGuess(double n, double T, const GasModule::GrainInterface&) const;

	/** Solves for the NLTE, given a total hydrogen density n, an initial (electron)
	    temperature guess, and a Spectrum object containing the radiation field in specific
	    intensity per frequency units. */
	Solution solveTemperature(double n, double Tinit, const Spectrum& specificIntensity,
	                          const GasModule::GrainInterface&) const;

	/** Calculate several extra contributions to the heating of the grains (collisions
	    (Draine and Bertoldi 1996), H2 formation on the surface (Takahashi 2001). Passing
	    some intermediary results of either the H2 level calculation or the grain
	    photoelectric effect calculation might help in speeding up this computation, as well
	    as using caching for the radiation emitted by the grains. */
	void updateGrainTemps(const Solution& s, const GasModule::GrainInterface&) const;

	/** Distills the solution object into the necessary information to retrieve opacity and
	    emissivity. Grids on which the opacity and emissivity will be discretized (before
	    being stored in the gas state) need to be provided. */
	GasModule::GasState makeGasStateFromSolution(const Solution&,
	                                             const Array& oFrequencyv,
	                                             const Array& eFrequencyv) const;

	/** Copies some values from the Solution, or recalculates them, and puts these in the
	    GasDiagnostics object */
	void fillGasDiagnosticsFromSolution(const Solution&, const GasModule::GrainInterface&,
	                                    GasDiagnostics*) const;

	/** Calculates all the densities for a fixed temperature. Is repeatedly called by this
	    class until equilibrium is found. */
	Solution solveDensities(double n, double T, const Spectrum& specificIntensity,
	                        const GasModule::GrainInterface&,
	                        const GasInterfaceImpl::Solution* previous = nullptr,
	                        double h2FormationOverride = -1) const;

	/** @name Properties of final state

	    These function calculate the properties which depend on the state of the gas. They
	    all take a \c Solution object as an argument, which can be obtained by calling \c
	    calculateDensities(). */
	/**@{*/
	/** The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emissivityv(const Solution&, const Array& eFrequencyv) const;

	/** The total opacity at each frequency in 1 / cm */
	Array opacityv(const Solution&, const Array& oFrequencyv) const;

	/** A possible scattering opacity */
	Array scatteringOpacityv(const Solution&, const Array& oFrequencyv) const;

	Array radiativeRecombinationEmissivityv(const Solution&,
	                                        const Array& eFrequencyv) const;

	Array freeFreeEmissivityv(const Solution&, const Array& eFrequencyv) const;

	Array lineEmissivityv(const Solution&, const Array& eFrequencyv) const;

	Array ionizationOpacityv(const Solution&, const Array& eFrequencyv) const;

	/** Total cooling */
	double cooling(const Solution&) const;

	/** The total heating, including the grain photoelectric effect, in erg / s / cm^3. */
	double heating(const Solution&, const GasModule::GrainInterface&) const;

	/** The heating by everything but the grains (relatively cheap) */
	double heating(const Solution&) const;

	/** The heating by the grains only (expensive to calculate) */
	double grainHeating(const Solution&, const GasModule::GrainInterface&) const;
	/**@}*/

	/** Easy access for some frequently used quantities */
	inline double nH(const Solution& s) const { return s.speciesNv[_inH]; }
	inline double nH2(const Solution& s) const { return s.speciesNv[_inH2]; }
	inline double np(const Solution& s) const { return s.speciesNv[_inp]; }
	inline double ne(const Solution& s) const { return s.speciesNv[_ine]; }

	/** Expression for the ionized fraction. */
	inline double f(const Solution& s) const
	{
		return s.speciesNv(_inp) / (s.speciesNv(_inH) + 2 * s.speciesNv(_inH2));
	}
	

private:
	// These are shorthand for ChemicalNetwork::speciesIndex.at["name"]
	int _ine, _inp, _inH, _inH2;
	std::unique_ptr<ChemistrySolver> _chemSolver;

	/* Pointers to other parts of the implementation, to make late initialization
	   possible */
	std::unique_ptr<HydrogenLevels> _atomicLevels;
	std::unique_ptr<H2Levels> _molecular;
	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound;
	std::unique_ptr<FreeFree> _freeFree;
};

#endif /* _HYDROGENCALCULATOR_H_ */
