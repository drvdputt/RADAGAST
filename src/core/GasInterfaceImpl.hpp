#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.hpp"
#include "GasSolution.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "GrainPhotoelectricEffect.hpp"
#include "NLevel.hpp"
#include "Spectrum.hpp"

class ChemistrySolver;
class FreeBound;
class FreeFree;
class GasDiagnostics;
class HydrogenLevels;
class H2Levels;
class SpeciesModelManager;

/** This class is the backbone of the whole code. It calculates the abundances, level
    populations and net heating rate of hydrogen repeatedly, until the equilibrium temperature
    is found.

    The parts of the code are grouped by the kind of process. When an instance of
    GasInterfaceImpl is constructed, a few members which provide implementations of the various
    processes are constructed too. Objects of the classes FreeFree, FreeBound and Chemistry are
    created, which decribe respectively the Brehmsstrahlung, the recombination continuum and the
    chemical network. During their construction, they each read the data they need to perform
    their functions.

    Another important member is the SpeciesModelManager, which contains all constant data for
    specialized species models (currently H and H2). I also allows this class to create HModel
    and H2Model objects, which will be given to a GasSolution object to store all the details of
    H and H2.

    The main methods of this class return a GasSolution object, which is given a callback
    pointer. Any quantities which depend on the specifics of the final equilibrium will be
    calculated in that class, with some callbacks to this class to access some functions which
    depend on constants stored here. */
class GasInterfaceImpl
{
public:
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
	GasSolution solveInitialGuess(double n, double T,
	                              const GasModule::GrainInterface&) const;

	/** Solves for the NLTE, given a total hydrogen density n, an initial (electron)
	    temperature guess, and a Spectrum object containing the radiation field in specific
	    intensity per frequency units. */
	GasSolution solveTemperature(double n, double Tinit, const Spectrum& specificIntensity,
	                             const GasModule::GrainInterface&) const;

	/** Calculates all the densities for a fixed temperature. Is repeatedly called by this
	    class until equilibrium is found. */
	GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity,
	                           const GasModule::GrainInterface&,
	                           const GasSolution* previous = nullptr,
	                           double h2FormationOverride = -1) const;

	/** Emission coefficient for the free-bound recombination continuum (per cm-6, need to
	    multiply with ne * np / 4pi) */
	Array radiativeRecombinationEmissivityv(double T, const Array& eFrequencyv) const;

	/** Emission coefficient for the free-free continuum (per cm-6, need to multiply with
	    ne*np / 4pi to obtain [erg s-1 cm-2 sr-1 Hz-1] */
	Array freeFreeEmissivityv(double T, const Array& eFrequencyv) const;

	/** Opacity coefficient for the free-free continuum (per cm-6, need to multiply with
	    ne*np to obtain [cm-1]) TODO: find out if this is relevant at all (I think it only
	    matters in radio, so maybe for 21 cm). */
	Array freeFreeOpacityv(double T, const Array& oFrequencyv) const;

	/** Cooling by free-free emission [erg s-1 cm-3] */
	double freeFreeCool(double np_ne, double T) const;

	/** Leaving this here for now, instead of moving it to GasSolution. GasSolution is not
	    the owner of the GrainInterface object, and it would be weird to modify it from
	    there. */
	void updateGrainTemps(const GasSolution& s,
			      GasModule::GrainInterface& g) const;


private:
	// These are shorthand for ChemicalNetwork::speciesIndex.at["name"]
	int _ine, _inp, _inH, _inH2;
	std::unique_ptr<ChemistrySolver> _chemSolver;

	/* Pointers to other parts of the implementation, to make late initialization
	   possible */
	std::unique_ptr<SpeciesModelManager> _manager;

	/* Continuum processes */
	std::unique_ptr<FreeBound> _freeBound;
	std::unique_ptr<FreeFree> _freeFree;
};

#endif /* _HYDROGENCALCULATOR_H_ */
