#include "GasInterface.hpp"
#include "GasInterfaceImpl.hpp"
#include <gsl/gsl_errno.h>
#include <algorithm>

namespace RADAGAST
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

    void GasInterface::updateGasState(GasState& gs, double n, const std::valarray<double>& meanIntensityv,
                                      double fshield, GrainInterface& gri, GasDiagnostics* gd) const
    {
        _pimpl->updateGasState(gs, n, meanIntensityv, fshield, gri, gd);
    }

    std::valarray<double> GasInterface::serialize(const GasState& gs) const
    {
        std::valarray<double> blob(6);
        blob[0] = gs._t;
        std::copy(std::begin(gs._nv), std::end(gs._nv), std::begin(blob) + 1);
        blob[blob.size() - 1] = gs._n2s;
        return blob;
    }

    void GasInterface::deserialize(GasState& gs, const std::valarray<double>& blob) const
    {
        // this is all just hardcoded. Assumes that blob and current size of _nv in gas state
        // are compatible. If not, will (and should) probably crash.
        gs.setMembers(blob[0], std::valarray<double>(&(blob[1]), gs._nv.size()), blob[blob.size() - 1]);
    }

    std::vector<std::string> GasInterface::serializationInfo() const
    {
        // TODO this should be made smarter in case the chemical network becomes configurable
        return {"T", "ne", "np", "nH", "nH2", "n2s"};
    }

    std::valarray<double> GasInterface::emissivityBasic(const GasState& gs, bool SI) const
    {
        return _pimpl->emissivityBasic(gs, SI);
    }

    std::valarray<double> GasInterface::opacityBasic(const GasState& gs, bool SI) const
    {
        return _pimpl->opacityBasic(gs, SI);
    }

    std::valarray<double> GasInterface::emissivityWithLines(const RADAGAST::GasState& gs,
                                                            const std::valarray<double>& meanIntensityv,
                                                            const GrainInterface& gri, bool SI, bool addHLines,
                                                            bool addH2Lines) const
    {
        return _pimpl->emissivityWithLines(gs, meanIntensityv, gri, SI, addHLines, addH2Lines);
    }

    std::valarray<double> GasInterface::opacityWithLines(const RADAGAST::GasState& gs,
                                                         const std::valarray<double>& meanIntensityv,
                                                         const GrainInterface& gri, bool SI, bool addHLines,
                                                         bool addH2Lines) const
    {
        return _pimpl->opacityWithLines(gs, meanIntensityv, gri, SI, addHLines, addH2Lines);
    }

    std::string GasInterface::quickInfo(const GasState& gs, const std::valarray<double>& meanIntensity) const
    {
        return _pimpl->quickInfo(gs, meanIntensity);
    }

    int GasInterface::index(const std::string& name) const { return _pimpl->index(name); }

    GasSolution GasInterface::solveTemperature(double n, const Spectrum& meanIntensity, double fshield,
                                               RADAGAST::GrainInterface& gri) const
    {
        return _pimpl->solveTemperature(n, meanIntensity, fshield, gri);
    }

    GasSolution GasInterface::solveDensities(double n, double T, const Spectrum& meanIntensity, double fshield,
                                             RADAGAST::GrainInterface& gri, double h2FormationOverride) const
    {
        return _pimpl->solveDensities(n, T, meanIntensity, fshield, gri, h2FormationOverride);
    }

    double GasInterface::solveDensities(GasSolution& s, double n, double T, const Spectrum& meanIntensity,
                                        bool startFromCurrent, double h2FormationOverride) const
    {
        return _pimpl->solveDensities(s, n, T, meanIntensity, startFromCurrent, h2FormationOverride);
    }

    // GasSolution GasInterface::solveDensitiesNoH2(double n, double T,
    //                                              const Spectrum& meanIntensity,
    //                                              RADAGAST::GrainInterface& gri) const
    // {
    // 	return _pimpl->solveDensitiesNoH2(n, T, meanIntensity, gri);
    // }

}  // namespace RADAGAST
