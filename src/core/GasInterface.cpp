#include "GasInterface.hpp"
#include "GasInterfaceImpl.hpp"

using namespace std;

namespace GasModule
{
GasInterface::GasInterface(const valarray<double>& iFrequencyv,
                           const valarray<double>& oFrequencyv,
                           const valarray<double>& eFrequencyv, const string& atomChoice,
                           const string& moleculeChoice)
                : _pimpl{make_unique<GasInterfaceImpl>(iFrequencyv, oFrequencyv, eFrequencyv,
                                                       atomChoice, moleculeChoice)}
{
}

GasInterface::~GasInterface() = default;

GasInterface::GasInterface(GasInterface&&) = default;

std::valarray<double> GasInterface::iFrequencyv() const { return _pimpl->iFrequencyv(); }

std::valarray<double> GasInterface::oFrequencyv() const { return _pimpl->oFrequencyv(); }

std::valarray<double> GasInterface::eFrequencyv() const { return _pimpl->eFrequencyv(); }

void GasInterface::updateGasState(GasState& gs, double n,
                                  const std::valarray<double>& specificIntensityv,
                                  GrainInterface& gri, GasDiagnostics* gd) const
{
	_pimpl->updateGasState(gs, n, specificIntensityv, gri, gd);
}

void GasInterface::initializeGasState(GasState& gs, double n, double T, GrainInterface& gri,
                                      GasDiagnostics* gd) const
{
	_pimpl->initializeGasState(gs, n, T, gri, gd);
}

double GasInterface::emissivity_SI(const GasState& gs, size_t iFreq) const
{
	return _pimpl->emissivity_SI(gs, iFreq);
}

double GasInterface::opacity_SI(const GasState& gs, size_t iFreq) const
{
	return _pimpl->opacity_SI(gs, iFreq);
}

GasSolution GasInterface::solveInitialGuess(double n, double T, GrainInterface& gri) const
{
	return _pimpl->solveInitialGuess(n, T, gri);
}

GasSolution GasInterface::solveTemperature(double n, const Spectrum& specificIntensity,
                                           GasModule::GrainInterface& gri) const
{
	return _pimpl->solveTemperature(n, specificIntensity, gri);
}

GasSolution GasInterface::solveDensities(double n, double T, const Spectrum& specificIntensity,
                                         GasModule::GrainInterface& gri,
                                         double h2FormationOverride) const
{
	return _pimpl->solveDensities(n, T, specificIntensity, gri, h2FormationOverride);
}

void GasInterface::solveDensities(GasSolution& s, double n, double T,
                                  const Spectrum& specificIntensity,
                                  GasModule::GrainInterface& gri, bool startFromCurrent,
                                  double h2FormationOverride) const
{
	_pimpl->solveDensities(s, n, T, specificIntensity, gri, startFromCurrent,
	                       h2FormationOverride);
}

// GasSolution GasInterface::solveDensitiesNoH2(double n, double T,
//                                              const Spectrum& specificIntensity,
//                                              GasModule::GrainInterface& gri) const
// {
// 	return _pimpl->solveDensitiesNoH2(n, T, specificIntensity, gri);
// }

} // namespace GasModule
