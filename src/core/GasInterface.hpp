#ifndef CORE_GASINTERFACE_HPP
#define CORE_GASINTERFACE_HPP

#include "GasState.hpp"
#include "GrainInterface.hpp"
#include <memory>
#include <string>

class GasDiagnostics;
class GasInterfaceImpl;
class GasSolution;
class Spectrum;

namespace GasModule
{
    /** The interface class that other codes should use. We use the PIMPL pattern here to minimize
        the amount of includes necessary on the client side, so that the internals of the gas
        module can be changed without having to recompile the client code. For the documentation,
        see GasInterfaceImpl. */
    class GasInterface
    {
    public:
        GasInterface(const std::valarray<double>& iFrequencyv, const std::valarray<double>& oFrequencyv,
                     const std::valarray<double>& eFrequencyv, const std::string& atomChoice = "",
                     const std::string& moleculeChoice = "");

        /** Needed because of the unique_ptr member. Also, putting "= default" here actually works
            with GDB if I remember correctly, but not with clang. In the cpp file, it works for
            both because there GasInterfaceImpl is a complete type. */
        ~GasInterface();

        /** This move constructor enables move semantics when utility functions in Testing return a
            GasInterface object. Again, necessary because of the unique_ptr. */
        GasInterface(GasInterface&&);

        const std::valarray<double>& iFrequencyv() const;
        const std::valarray<double>& oFrequencyv() const;
        const std::valarray<double>& eFrequencyv() const;

        void updateGasState(GasState&, double n, const std::valarray<double>& specificIntensityv, GrainInterface& gri,
                            GasDiagnostics* gd = nullptr) const;
        void initializeGasState(GasState&, double n, double T, GrainInterface&, GasDiagnostics* gd = nullptr) const;
        double emissivity_SI(const GasState& gs, size_t iFreq) const;
        double opacity_SI(const GasState& gs, size_t iFreq) const;

        GasSolution solveInitialGuess(double n, double T, GrainInterface&) const;
        GasSolution solveTemperature(double n, const Spectrum& specificIntensity, GasModule::GrainInterface&) const;
        GasSolution solveDensities(double n, double T, const Spectrum& specificIntensity, GasModule::GrainInterface&,
                                   double h2FormationOverride = -1) const;
        void solveDensities(GasSolution&, double n, double T, const Spectrum& specificIntensity,
                            bool startFromCurrent = false, double h2FormationOverride = -1) const;
        GasSolution solveDensitiesNoH2(double n, double T, const Spectrum& specificIntensity,
                                       GasModule::GrainInterface&) const;

    private:
        std::unique_ptr<GasInterfaceImpl> _pimpl;
    };
}  // namespace GasModule
#endif  // CORE_GASINTERFACE_HPP
