#include "GasInterface.hpp"
#include "GasInterfaceImpl.hpp"
#include <gsl/gsl_errno.h>

namespace GasModule
{
    namespace
    {
        gsl_error_handler_t* _oldGSLErrorHandler;
    }
    void GasInterface::errorHandlersOff() { _oldGSLErrorHandler = gsl_set_error_handler_off(); }

    GasInterface::GasInterface(const std::valarray<double>& iFrequencyv, const std::valarray<double>& oFrequencyv,
                               const std::valarray<double>& eFrequencyv)
        : _pimpl{std::make_unique<GasInterfaceImpl>(iFrequencyv, oFrequencyv, eFrequencyv)}
    {}

    GasInterface::~GasInterface() = default;

    GasInterface::GasInterface(GasInterface&&) = default;

    const std::valarray<double>& GasInterface::iFrequencyv() const { return _pimpl->iFrequencyv(); }

    const std::valarray<double>& GasInterface::oFrequencyv() const { return _pimpl->oFrequencyv(); }

    const std::valarray<double>& GasInterface::eFrequencyv() const { return _pimpl->eFrequencyv(); }

    void GasInterface::updateGasState(GasState& gs, double n, const std::valarray<double>& specificIntensityv,
                                      GrainInterface& gri, GasDiagnostics* gd) const
    {
        _pimpl->updateGasState(gs, n, specificIntensityv, gri, gd);
    }

    std::valarray<double> GasInterface::emissivityBasic(const GasState& gs, bool SI) const
    {
        return _pimpl->emissivityBasic(gs, SI);
    }

    std::valarray<double> GasInterface::opacityBasic(const GasState& gs, bool SI) const
    {
        return _pimpl->opacityBasic(gs, SI);
    }

    std::valarray<double> GasInterface::emissivityWithLines(const GasModule::GasState& gs,
                                                            const std::valarray<double>& specificIntensityv,
                                                            const GrainInterface& gri, bool SI, bool addHLines,
                                                            bool addH2Lines) const
    {
        return _pimpl->emissivityWithLines(gs, specificIntensityv, gri, SI, addHLines, addH2Lines);
    }

    std::valarray<double> GasInterface::opacityWithLines(const GasModule::GasState& gs,
                                                         const std::valarray<double>& specificIntensityv,
                                                         const GrainInterface& gri, bool SI, bool addHLines,
                                                         bool addH2Lines) const
    {
        return _pimpl->opacityWithLines(gs, specificIntensityv, gri, SI, addHLines, addH2Lines);
    }

    std::string GasInterface::quickInfo(const GasState& gs, const std::valarray<double>& specificIntensity) const
    {
        return _pimpl->quickInfo(gs, specificIntensity);
    }

    int GasInterface::index(const std::string& name) const { return _pimpl->index(name); }

    GasSolution GasInterface::solveTemperature(double n, const Spectrum& specificIntensity,
                                               GasModule::GrainInterface& gri) const
    {
        return _pimpl->solveTemperature(n, specificIntensity, gri);
    }

    GasSolution GasInterface::solveDensities(double n, double T, const Spectrum& specificIntensity,
                                             GasModule::GrainInterface& gri, double h2FormationOverride) const
    {
        return _pimpl->solveDensities(n, T, specificIntensity, gri, h2FormationOverride);
    }

    double GasInterface::solveDensities(GasSolution& s, double n, double T, const Spectrum& specificIntensity,
                                        bool startFromCurrent, double h2FormationOverride) const
    {
        return _pimpl->solveDensities(s, n, T, specificIntensity, startFromCurrent, h2FormationOverride);
    }

    // GasSolution GasInterface::solveDensitiesNoH2(double n, double T,
    //                                              const Spectrum& specificIntensity,
    //                                              GasModule::GrainInterface& gri) const
    // {
    // 	return _pimpl->solveDensitiesNoH2(n, T, specificIntensity, gri);
    // }

}  // namespace GasModule
