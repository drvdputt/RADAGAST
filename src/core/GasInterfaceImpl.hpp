#ifndef CORE_GASINTERFACEIMPL_HPP
#define CORE_GASINTERFACEIMPL_HPP

#include "Array.hpp"
#include "FreeBound.hpp"
#include "FreeFree.hpp"
#include "GasSolution.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "SimpleHChemistry.hpp"
#include "SpeciesModelManager.hpp"
#include "Spectrum.hpp"

class Chemistry;

/** This class is the backbone of the whole code. It calculates the abundances, level populations
    and net heating rate of hydrogen repeatedly, until the equilibrium temperature is found.

    When an instance is constructed, a few members which provide implementations of the various
    processes are constructed too. FreeFree, FreeBound and Chemistry objects are created, which
    decribe respectively the Brehmsstrahlung, the recombination continuum and the chemical network.
    During their construction, they each read the data they need to perform their functions.

    The main methods of this class return a GasSolution object. At construction, the GasSolution
    takes ownership of HModel and H2Model objects created using the SpeciesModelManager, and
    receives (const) references to the FreeFree and FreeBound instances contained here. From the
    GasSolution, full details about the solution can be obtained (see its documentation).

    From a GasSolution object, more lightweight GasState objects can also be obtained. A few
    convenience functions are provided in this class to skip having to make a GasSolution object (an
    intermediary GasSolution still appears though). Only using GasState objects drastically reduces
    the amount of includes needed. The GasState objects can be given to a few functions in this
    class to calculate the emissivity and opacity.

    All the functions in this class are reentrant, which means that they can safely be called
    simultaneously by different threads. All the non-constant variables are contained in the
    GasSolution objects, of which multiple can exist at the same time. */
class GasInterfaceImpl
{
public:
    /** Creates an instance of the gas module. Multiple frequency grids are used, which need to be
        specified by the user. Some configuration options in the form of strings are also provided.
        Currently, they only influence the settings of the H and H2 models, see
        SpeciesModelManager. */
    GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv, const Array& eFrequencyv,
                     const std::string& atomChoice = "", const std::string& moleculeChoice = "");

    /** The grid used to discretize the input radiation field */
    const std::valarray<double>& iFrequencyv() const { return _iFrequencyv; }

    /** The grid that will be used to discretize the output opacity. This is typically coarser
        because a radiative transfer algorithm usually needs the opacity in each grid cell. */
    const std::valarray<double>& oFrequencyv() const { return _oFrequencyv; }

    /** The grid on which the emissivity is calculated. */
    const std::valarray<double>& eFrequencyv() const { return _eFrequencyv; }

    /** The most convenient way to run the code for a cell. A minimal set of results is stored in
        the given GasState object. The exact contents of the GasState are not known to the user.
        Through this interface, the variable information contained in the gas state can be combined
        with other constants and functions to retrieve the opacity and emissivity of the gas.

        The other arguments reflect the physical conditions in the cell for which a gas state is
        being calculated. The density of hydrogen nuclei, n; the ambient radiation field in [erg s-1
        cm-1 Hz-1] units, as discretized on the @c iFrequencyv grid given as an argument to the
        constructor, originally; and a GrainInterface object, describing the properties of the
        grains within the cell for which the user wants to solve the gas state. See the @c
        GrainInterface documentation for information on how to construct one of these objects. The
        absorption efficiency of the grains currently needs to be tabulated for the frequencies
        contained in iFrequencyv.

        Note that the GrainInterface instance can be modified; the temperatures are recalculated to
        take into account the effect of gas-grain collisions. */
    void updateGasState(GasModule::GasState&, double n, const std::valarray<double>& specificIntensityv,
                        GasModule::GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

    /** Does the same as the above, but without an input radiation field. Instead, a blackbody of
        the given temperature is used to calculate the GasState. It is recommended to apply this
        function to all gas states before starting a simulation. */
    void initializeGasState(GasModule::GasState&, double n, double T, GasModule::GrainInterface&,
                            GasDiagnostics* gd = nullptr) const;

    /** The emissivity in SI units, for a given frequency index. (converted using 1 erg cm-3 s-1
        Hz-1 sr-1 = 0.1 J m-3 s-1 Hz-1 sr-1). [W m-3 Hz-1 sr-1] */
    double emissivity_SI(const GasModule::GasState& gs, size_t iFreq) const;

    /** The total opacity in SI units (converted from cm-1 = 100 * m-1). [m-1] */
    double opacity_SI(const GasModule::GasState& gs, size_t iFreq) const;

    /** This convenience function runs solveTemperature for a blackbody radiation field. The final
        temperature is usually close to the given color temperature T. */
    GasSolution solveInitialGuess(double n, double T, GasModule::GrainInterface&) const;

    /** Find a fully self-consistent solution and temperature, given a total hydrogen density n and
        a Spectrum object describing the radiation field in specific intensity units [erg s-1 sr-1
        cm-2]. */
    GasSolution solveTemperature(double n, const Spectrum& specificIntensity, GasModule::GrainInterface&) const;

    /** Find a fully self-consistent solution for a fixed temperature (heating/cooling will be out
        of equilibrium). */
    GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity, GasModule::GrainInterface&,
                               double h2FormationOverride = -1) const;

    /** Recalculate the densities for a GasSolution object. If startFromCurrent is true, the
        current contents of the GasSolution are used as an initial guess if the given temperature
        doesn't differ too much (< factor 2) from the previous temperature. Otherwise, an initial
        guess is made for the chemistry, where the ionized fraction is based on the radiation
        field, and the initial molecular fraction is 0.1. */
    void solveDensities(GasSolution&, double n, double T, const Spectrum& specificIntensity,
                        bool startFromCurrent = false, double h2FormationOverride = -1) const;

    /** NOT IMPLEMENTED. Solve for a fixed temperature, forcing H2 to zero. Useful for
        ionization-only tests. */
    GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& specificIntensity,
                                   GasModule::GrainInterface&) const;

private:
    /** Construct a new GasSolution model, passing all the necessary (references to) objects. The
        SpeciesModelManager is used to create the HModel and H2Model. */
    GasSolution makeGasSolution(const Spectrum& specificIntensity, const GasModule::GrainInterface*) const;

    /** Create a species vector which is zero everywhere, except for e-, p+, H, and H2. Arguments:
        the total amount of H nuclei n, the ionized to total fraction (np / n), the molecular to
        neutral fraction (2 * nH2 / (nH + 2 * nH2)). */
    EVector guessSpeciesNv(double n, double ionToTotalFrac, double moleculeToNeutralFrac) const;

    std::valarray<double> _iFrequencyv;
    std::valarray<double> _oFrequencyv;
    std::valarray<double> _eFrequencyv;

    SimpleHChemistry _chemistry{};
    SpeciesModelManager _manager;
    FreeBound _freeBound{};
    FreeFree _freeFree{};
};

#endif  // CORE_GASINTERFACEIMPL_HPP
