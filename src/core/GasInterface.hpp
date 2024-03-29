#ifndef CORE_GASINTERFACE_HPP
#define CORE_GASINTERFACE_HPP

#include "GasState.hpp"
#include "GrainInterface.hpp"
#include <memory>
#include <string>

namespace RADAGAST
{
    class GasDiagnostics;
    class GasInterfaceImpl;
    class GasSolution;
    class Spectrum;

    /** This is the interface class that other codes should use. We use the PIMPL pattern here
        to minimize the amount of includes necessary on the client side, so that the internals
        of the gas module can be changed without having to recompile large parts of the client
        code. For implementation details, GasInterfaceImpl. Another advantage is that when
        installing RADAGAST, only a handful of headers need to be copied to the include/
        directory. Including GasInterfaceImpl.hpp in the client code directly, would pull in a
        lot of extra headers.

        There are several groups of functions available in this interface. First, there's the
        constructor and several getters for the members that were set during construction. Any data
        given during construction are considered constant during the simulation, and some
        precalculations might be performed based on them.

        The second group consists of functions for working with @c GasState objects. These
        provide a lightweight way to store the result of the calculation, but they are agnostic
        of any details about the model. @c updateGasState and @c initializeGasState change the
        gas state. The rest extract information from the gas state, either directly, or by using
        some constants stored in the GasInterface instance, or by performing some
        (re)calculations.

        The third group of functions gives full access to the functions that create and modify @c
        GasSolution objects. These are heavyweight, containing all the details of the equilibrium
        calculation. Using them requires including the @c GasSolution header, which in turn pulls
        in a lot of other classes. But long as the user does not need GasSolution objects, these
        functions can be ignored and the necessary includes are kept to a minimum: @c GasInterface,
        @c GrainInterface and @c GasState. */
    class GasInterface
    {
    public:
        // SETUP //

        /** Turn off any error handlers which could halt or abort the program. Currently this
            simply turns off the GSL error handler. If this code would become dependent on other
            libraries that have their own error handlers, this function would be a useful way to
            turn them off too. This static function has to be called before any other functionality
            of the gas module is used. */
        static void errorHandlersOff();

        /** Create an instance of the gas module. Multiple frequency grids are used, which need
            to be specified by the user. */
        GasInterface(const std::valarray<double>& iFrequencyv, const std::valarray<double>& oFrequencyv,
                     const std::valarray<double>& eFrequencyv);

        /** Defer the implementation of the destructor to the cpp file. Needed because of the
            unique_ptr member. Also, putting "= default" here actually works with GDB if I remember
            correctly, but not with clang. In the cpp file, it works for both because there
            GasInterfaceImpl is a complete type. */
        ~GasInterface();

        /** Enable move semantics when utility functions in Testing return a GasInterface object.
            Again, necessary because of the unique_ptr. */
        GasInterface(GasInterface&&);

        /** Return the frequency grid for the input radiation field. */
        const std::valarray<double>& iFrequencyv() const;

        /** Return the frequency grid for the output opacity. */
        const std::valarray<double>& oFrequencyv() const;

        /** Return the frequency grid for the output emissivity. */
        const std::valarray<double>& eFrequencyv() const;

        // UPDATE GASSTATE AND RETRIEVE INFORMATION //

        /** Update the contents of the given @c GasState object with the results of a new
            equilibrium calculation. The exact contents of the GasState are not known to the user.
            The information contained in the gas state can be combined with other constants and
            functions to retrieve the opacity and emissivity of the gas, as well as other
            diagnostices.

            The arguments specify the physical conditions in the cell for which the equilibrium is
            being calculated. The calculation requires: the density of hydrogen nuclei, n; the
            ambient radiation field in [erg s-1 cm-1 Hz-1] units, as discretized on the @c
            iFrequencyv grid given as an argument to the constructor, originally; and a @c
            GrainInterface object, describing the properties of the grains within the cell. See the
            @c GrainInterface documentation for information on how to construct one of these
            objects. The absorption efficiency of the grains currently needs to be tabulated for
            the same frequencies contained in iFrequencyv.

            Since H2 self-shielding is a non-local effect, it cannot be taken into account by
            RADAGAST (except with extremely high spectral resolution). The fshield argument
            provides a way to let the client provide an appropriate self-shielding factor
            (between 0 and 1), typically derived from its own radiative transfer calculations or
            approximation for NH2. If this factor is provided, the 'naive' free-space
            dissociation rate will be multiplied by it. (D = f * average FUV * integrated
            sigma).

            The last argument is an optional pointer to a GasDiagnostics object. If provided, the
            object will be filled with various diagnostic data (slow!). I'm not sure if I'm going
            to keep this functionality, now that the functions for working with @c GasSolution
            exist. */
        void updateGasState(GasState&, double n, const std::valarray<double>& meanIntensityv, double fshield,
                            GrainInterface& gri, GasDiagnostics* gd = nullptr) const;

        /** Serialize GasState into array of doubles, which makes it easier for the client code
            to implement storage and e.g. MPI communication. Not implemented in GasState,
            because GasState is agnostic about the details of the model, which we might need. */
        std::valarray<double> serialize(const GasState&, bool SI = false) const;

        /** Write serialized data into GasState. Typical use case:
            1. run RADAGAST in parallel (MPI)
            2. store gas states as plain doubles using serialize
            3. synchronize gas state data between processes using simple MPI calls
            4. deserialize the blobs and pass the resulting GasStates to RADAGAST to access gas information */
        void deserialize(GasState&, const std::valarray<double>& blob, bool SI = false) const;

        /** Describes every element of the output of serialize, in the form of one string per
            element. The output looks like {"T", "ne", "np", "nH", "nH2", others...}. GasState
            is agnostic of the density vector, but GasInterface has access to the chemical
            model, which is why this is implemented here. */
        std::vector<std::string> serializationInfo() const;

        /** Return the emissivity of the gas, discretized on @c eFrequencyv [erg cm-3 s-1 Hz-1
            sr-1]. If the second argument is @c true, SI units are used instead [W m-3 Hz-1
            sr-1]. */
        std::valarray<double> emissivityBasic(const GasState& gs, bool SI = false) const;

        /** Return the opacity of the gas [cm-1], discretized on @c oFrequencyv. If the second
            argument is @c true, SI units are used instead [m-1]. */
        std::valarray<double> opacityBasic(const GasState& gs, bool SI = false) const;

        /** Same as the above, but will add the emissivity of the lines. Some recalculations are
            necessary, hence the radiation field and grain details are needed again. */
        std::valarray<double> emissivityWithLines(const RADAGAST::GasState& gs,
                                                  const std::valarray<double>& meanIntensityv,
                                                  const GrainInterface& gri, bool SI, bool addHLines,
                                                  bool addH2Lines) const;

        /** Same as the above, but will add the opacity of the lines. Some recalculations are
            necessary, hence the radiation field and grain details are needed again. */
        std::valarray<double> opacityWithLines(const RADAGAST::GasState& gs,
                                               const std::valarray<double>& meanIntensityv, const GrainInterface& gri,
                                               bool SI, bool addHLines, bool addH2Lines) const;

        /** Create a string containing a 1-line overview of the result. Units are included in the
            string. */
        std::string quickInfo(const GasState& gs, const std::valarray<double>& meanIntensity) const;

        /** Return the index of the stored density in the gas state for a given species
            (specified as a string). If the given species is not available, -1 will be returned.
            The @c GasState itself has no knowledge of what its stored densities mean. The
            latter can be different depending on how the chemical network was set up. This index
            can be used to retrieve the desired density from a @c GasState object that was
            created through this @c GasInterface instance. */
        int index(const std::string& name) const;

        // ACCESS TO GASSOLUTION //

        /** Find the equilibrium temperature, by repeatedly solving the densities and evaluating th
            heating and cooling. @c updateGasState works via this function under the hood. */
        GasSolution solveTemperature(double n, const Spectrum& meanIntensity, double fshield,
                                     RADAGAST::GrainInterface&) const;

        /** Find the equilibrium densities for a fixed temperature (heating/cooling will be out of
            equilibrium). */
        GasSolution solveDensities(double n, double T, const Spectrum& meanIntensity, double fshield,
                                   RADAGAST::GrainInterface&, double h2FormationOverride = -1) const;

        /** Recalculate the densities for an existing GasSolution object. If startFromCurrent is
            true, the current contents of the GasSolution are potentially used as an initial guess
            If the given temperature is very different (> factor 2) from the temperature currently
            contained in the @c GasSolution, a more heuristic initial guess is made for the
            chemistry, where the ionized fraction is based on the radiation field, and the
            molecular fraction is 0.1. Returns the net heating (heating - cooling).*/
        double solveDensities(GasSolution&, double n, double T, const Spectrum& meanIntensity,
                              bool startFromCurrent = false, double h2FormationOverride = -1) const;

        /** NOT IMPLEMENTED. Solve for a fixed temperature, forcing H2 to zero. Useful for
            ionization-only tests. */
        GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& meanIntensity,
                                       RADAGAST::GrainInterface&) const;

    private:
        std::unique_ptr<GasInterfaceImpl> _pimpl;
    };
}
#endif  // CORE_GASINTERFACE_HPP
