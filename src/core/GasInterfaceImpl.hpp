#ifndef CORE_GASINTERFACEIMPL_HPP
#define CORE_GASINTERFACEIMPL_HPP

#include "Array.hpp"
#include "ChemistrySolver.hpp"
#include "FreeBound.hpp"
#include "FreeFree.hpp"
#include "GasSolution.hpp"
#include "GasState.hpp"
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
	/** Creates an instance of the gas module. Multiple frequency grids are used, which need
	    to be specified by the user. Some configuration options in the form of strings are
	    also provided. Currently, they only influence the settings of the H and H2
	    models. */
	GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv,
	                 const Array& eFrequencyv, const std::string& atomChoice = "",
	                 const std::string& moleculeChoice = "");

	/** The grid used to discretize the input radiation field */
	std::valarray<double> iFrequencyv() const { return _iFrequencyv; }

	/** The grid that will be used to discretize the output opacity. This is typically
	    coarser because a radiative transfer algorithm usually needs the opacity in each
	    grid cell. */
	std::valarray<double> oFrequencyv() const { return _oFrequencyv; }

	/** The grid on which the emissivity is calculated. */
	std::valarray<double> eFrequencyv() const { return _eFrequencyv; }

	// GASSTATE TOOLS (minimal information, works with GasSolution under the hood) //

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
	void updateGasState(GasModule::GasState&, double n,
	                    const std::valarray<double>& specificIntensityv,
	                    GasModule::GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

	/** Does the same as the above, but without an input radiation field. Instead, a
	    blackbody of the given temperature is used to calculate GasState. It is recommended
	    to apply this function to all gas states before starting a simulation. */
	void initializeGasState(GasModule::GasState&, double n, double T,
	                        GasModule::GrainInterface&, GasDiagnostics* gd = nullptr) const;

	/** The emissivity in SI units, for a given frequency index. (converted using 1 erg cm-3
	    s-1 Hz-1 sr-1 = 0.1 J m-3 s-1 Hz-1 sr-1). [W m-3 Hz-1 sr-1] */
	double emissivity_SI(const GasModule::GasState& gs, size_t iFreq) const;

	/** Returns the total opacity in SI units (converted from cm-1 = 100 * m-1). [m-1]*/
	double opacity_SI(const GasModule::GasState& gs, size_t iFreq) const;

	// GASSOLUTION TOOLS (full information) //

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

	std::valarray<double> _iFrequencyv;
	std::valarray<double> _oFrequencyv;
	std::valarray<double> _eFrequencyv;

	int _ine, _inp, _inH, _inH2;
	ChemistrySolver _chemSolver;
	SpeciesModelManager _manager;
	FreeBound _freeBound{};
	FreeFree _freeFree{};
};

#endif // CORE_GASINTERFACEIMPL_HPP
