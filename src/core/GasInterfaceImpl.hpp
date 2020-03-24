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

namespace GasModule
{
    class Chemistry;

    /** This class is the backbone of the whole code. It calculates the abundances, level
        populations and net heating rate of hydrogen repeatedly, until the equilibrium temperature
        is found. A public front-end to this class exists (\c GasInterface). The documentation for
        functions not documented here can be found there.

        When an instance is constructed, a few members which provide implementations of the various
        processes are constructed too. FreeFree, FreeBound and Chemistry objects are created, which
        decribe respectively the Brehmsstrahlung, the recombination continuum and the chemical
        network. During their construction, they each read the data they need to perform their
        functions.

        The main methods of this class return a GasSolution object. At construction, the
        GasSolution takes ownership of HModel and H2Model objects created using the
        SpeciesModelManager, and receives (const) references to the FreeFree and FreeBound
        instances contained here. From the GasSolution, full details about the solution can be
        obtained (see its documentation).

        From a GasSolution object, more lightweight GasState objects can also be obtained. A few
        convenience functions are provided in this class to skip having to make a GasSolution
        object (an intermediary GasSolution still appears though). Only using GasState objects
        drastically reduces the amount of includes needed. The GasState objects can be given to a
        few functions in this class to calculate the emissivity and opacity.

        All the functions in this class are reentrant, which means that they can safely be called
        simultaneously by different threads. All the non-constant variables are contained in the
        GasSolution objects, of which multiple can exist at the same time. */
    class GasInterfaceImpl
    {
    public:
        GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv, const Array& eFrequencyv,
                         const std::string& atomChoice = "", const std::string& moleculeChoice = "");

        const std::valarray<double>& iFrequencyv() const { return _iFrequencyv; }
        const std::valarray<double>& oFrequencyv() const { return _oFrequencyv; }
        const std::valarray<double>& eFrequencyv() const { return _eFrequencyv; }

        void updateGasState(GasModule::GasState&, double n, const std::valarray<double>& specificIntensityv,
                            GasModule::GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

        Array emissivity(const GasModule::GasState& gs, bool SI = false) const;
        Array opacity(const GasModule::GasState& gs, bool SI = false) const;
        std::string quickInfo(const GasModule::GasState& gs, const std::valarray<double>& specificIntensity) const;
        int index(const std::string& name) const;

        GasSolution solveTemperature(double n, const Spectrum& specificIntensity, GasModule::GrainInterface&) const;
        GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity, GasModule::GrainInterface&,
                                   double h2FormationOverride = -1) const;
        double solveDensities(GasSolution&, double n, double T, const Spectrum& specificIntensity,
                              bool startFromCurrent = false, double h2FormationOverride = -1) const;
        GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& specificIntensity,
                                       GasModule::GrainInterface&) const;

    private:
        /** Construct a new GasSolution model, passing all the necessary (references to) objects.
            The SpeciesModelManager is used to create the HModel and H2Model. */
        GasSolution makeGasSolution(const Spectrum& specificIntensity, const GasModule::GrainInterface*) const;

        /** Create a species vector which is zero everywhere, except for e-, p+, H, and H2.
            Arguments: the total amount of H nuclei n, the ionized to total fraction (np / n), the
            molecular to neutral fraction (2 * nH2 / (nH + 2 * nH2)). */
        EVector guessSpeciesNv(double n, double ionToTotalFrac, double moleculeToNeutralFrac) const;

        /** Utility function to get the product n_p n_e from a gas state. */
        double npne(const GasModule::GasState&) const;
        double nH(const GasModule::GasState&) const;

        /** Interpret the densities stored in the gas state as a species vector for the chemistry
            network contained in this class. */
        SpeciesVector speciesVector(const GasModule::GasState&) const;

        std::valarray<double> _iFrequencyv;
        std::valarray<double> _oFrequencyv;
        std::valarray<double> _eFrequencyv;

        SimpleHChemistry _chemistry{};
        SpeciesModelManager _manager;
        FreeBound _freeBound{};
        FreeFree _freeFree{};

        // Standalone, averaged (on ofrequencyv) H2 cross section. In the ideal case, the H2
        // cross section is calculated from the levels, but level-dependent stuff is not
        // available at the interface yet.
        Array _h2crossv;
    };
}
#endif  // CORE_GASINTERFACEIMPL_HPP
