#include "GasInterface.hpp"
#include "Array.hpp"
#include "Error.hpp"
#include "GasInterfaceImpl.hpp"
#include "Timer.hpp"

using namespace std;

namespace GasModule
{
GasInterface::GasInterface(const valarray<double>& iFrequencyv,
                           const valarray<double>& oFrequencyv,
                           const valarray<double>& eFrequencyv, const string& atomChoice,
                           const string& moleculeChoice)
                : _iFrequencyv{iFrequencyv}, _oFrequencyv{oFrequencyv},
                  _eFrequencyv{eFrequencyv}, _pimpl{make_unique<GasInterfaceImpl>(
                                                             atomChoice, moleculeChoice)}
{
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gasState, double n, double Tinit,
                                  const valarray<double>& specificIntensityv,
                                  GrainInterface& grainInfo, GasDiagnostics* gd) const
{
	Timer t("Update gas state");
	// Create a spectrum object which makes it easier to pass around the frequencies and the
	// values for the specific intensity. It can also be used to easily interpolate the
	// specific intensity for a certain frequency.
	Spectrum specificIntensity(_iFrequencyv, specificIntensityv);
	if (n > 0)
	{
		GasSolution s = _pimpl->solveTemperature(n, Tinit, specificIntensity, grainInfo);
		// To fix the temperature, use this:
		// GasInterfaceImpl::Solution s = _pimpl->solveDensities(n, 6000., specificIntensity, grainInfo);
		gasState = s.makeGasState(_oFrequencyv, _eFrequencyv);
		if (gd)
			s.fillDiagnostics(gd);
	}
	else
		zeroOpticalProperties(gasState);
}

void GasInterface::initializeGasState(GasState& gasState, double n, double T,
                                      GrainInterface& grainInfo, GasDiagnostics* gd) const
{
	if (n > 0)
	{
		GasSolution s = _pimpl->solveInitialGuess(n, T, grainInfo);
		gasState = s.makeGasState(_oFrequencyv, _eFrequencyv);
		if (gd)
			s.fillDiagnostics(gd);
	}
	else
		zeroOpticalProperties(gasState);
}

double GasInterface::emissivity_SI(const GasState& gs, size_t iFreq) const
{
	return 0.1 * gs._emissivityv[iFreq];
}

// 1 / cm = 100 / m
double GasInterface::opacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * (gs._opacityv[iFreq] + gs._scatteringOpacityv[iFreq]);
}

double GasInterface::scatteringOpacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * gs._scatteringOpacityv[iFreq];
}
double GasInterface::absorptionOpacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * gs._opacityv[iFreq];
}

void GasInterface::zeroOpticalProperties(GasState& gs) const
{
	gs._emissivityv = Array(_eFrequencyv.size());
	gs._opacityv = Array(_oFrequencyv.size());
	gs._scatteringOpacityv = Array(_oFrequencyv.size());
	gs._temperature = 0;
	gs._densityv = Array{0, 0, 0, 0};
}
} /* namespace GasModule */
