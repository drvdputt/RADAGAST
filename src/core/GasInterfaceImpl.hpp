#ifndef CORE_GASINTERFACEIMPL_HPP
#define CORE_GASINTERFACEIMPL_HPP

#include "Array.hpp"
#include "FreeBound.hpp"
#include "FreeFree.hpp"
#include "GasSolution.hpp"
#include "GrainInterface.hpp"
#include "SpeciesModelManager.hpp"
#include "Spectrum.hpp"

class ChemistrySolver;

/** This class is the backbone of the whole code. It calculates the abundances, level
    populations and net heating rate of hydrogen repeatedly, until the equilibrium temperature
    is found.

    When an instance is constructed, a few members which provide implementations of the various
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
    calculated in that class. This class has a handful of extra functions, to give the
    GasSolution access to some of the functionality of FreeFree and FreeBound. */
class GasInterfaceImpl
{
public:
	/** The constructor takes the same arguments as GasInterface. They are simply passed to
	    the SpeciesModelManager. */
	GasInterfaceImpl(const std::string& atomChoice, const std::string& moleculeChoice);

	~GasInterfaceImpl();

	/** This convenience function runs solveTemperature for a blackbody radiation field. The
	    final temperature should be close to the given blackbody temperature (?) */
	GasSolution solveInitialGuess(double n, double T, GasModule::GrainInterface&) const;

	/** Find a fully self-consistent solution, given a total hydrogen density n and a
	    Spectrum object describing the radiation field in specific intensity units [erg s-1
	    sr-1 cm-2]. */
	GasSolution solveTemperature(double n, const Spectrum& specificIntensity,
	                             GasModule::GrainInterface&) const;

	/** Calculates all the densities for a fixed temperature. */
	GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity,
	                           GasModule::GrainInterface&,
	                           double h2FormationOverride = -1) const;

	/** Recalculate the densities for a GasSolution object, with a different temperature. Is
	    repeatedly called by this class until equilibrium is found. */
	void solveDensities(GasSolution&, double n, double T, const Spectrum& specificIntensity,
	                    GasModule::GrainInterface&, bool startFromCurrent = false,
	                    double h2FormationOverride = -1) const;

	/** TODO: rework old code */
	GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& specificIntensity,
	                               GasModule::GrainInterface&) const;

	/** Calculate several extra contributions to the heating of the grains (collisions
	    (Draine and Bertoldi 1996), H2 formation on the surface (Takahashi 2001). Passing
	    some intermediary results of either the H2 level calculation or the grain
	    photoelectric effect calculation might help in speeding up this computation, as well
	    as using caching for the radiation emitted by the grains.

	    Leaving this here for now, instead of moving it to GasSolution. GasSolution is not
	    the owner of the GrainInterface object, and it would be weird to modify it from
	    there. */
	void updateGrainTemps(const GasSolution& s, GasModule::GrainInterface& g) const;

private:
	GasSolution makeGasSolution(const Spectrum& specificIntensity,
	                            const GasModule::GrainInterface&) const;
	EVector guessSpeciesNv(double n, double ionToTotalFrac,
	                       double moleculeToNeutralFrac) const;

	int _ine, _inp, _inH, _inH2;
	std::unique_ptr<ChemistrySolver> _chemSolver;

	/* Pointers to other parts of the implementation, to make late initialization
	   possible */
	SpeciesModelManager _manager;

	/* Continuum processes (no arguments needed) */
	FreeBound _freeBound{};
	FreeFree _freeFree{};
};

#endif // CORE_GASINTERFACEIMPL_HPP
