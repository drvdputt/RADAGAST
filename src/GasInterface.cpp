#include "GasInterface.h"
#include "Array.h"

#include "GasInterfaceImpl.h"

using namespace std;

GasInterface::GasInterface(const valarray<double>& frequencyv)
                : _frequencyv(frequencyv), _pimpl(make_unique<GasInterfaceImpl>(frequencyv))
{
}

GasInterface::GasInterface(const valarray<double>& frequencyv, bool improveGrid)
                : _frequencyv(frequencyv),
                  _pimpl(make_unique<GasInterfaceImpl>(frequencyv, improveGrid))
{
	_frequencyv = _pimpl->frequencyv();
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gs, double density, double Tinit,
                                  const valarray<double>& specificIntensityv) const
{
	gs._previousISRFv = specificIntensityv;
	if (density > 0)
		_pimpl->solveBalance(gs, density, Tinit, specificIntensityv);
	else
		zeroOpticalProperties(gs);
}

void GasInterface::initializeGasState(GasState& gs, double density, double temperature) const
{
	gs._previousISRFv = Array(_frequencyv.size());
	if (density > 0)
		_pimpl->solveInitialGuess(gs, density, temperature);
	else
		zeroOpticalProperties(gs);
}

double GasInterface::effectiveEmissivity_SI(const GasState& gs, size_t iFreq) const
{

	double r = 0.1 * (gs._emissivityv[iFreq] /*-
	                  gs._scatteringOpacityv[iFreq] * gs._previousISRFv[iFreq]*/);
	return r > 0 ? r : 0;
}

// 1 / cm = 100 / m
double GasInterface::opacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * gs._opacityv[iFreq];
}
double GasInterface::scatteringOpacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * gs._scatteringOpacityv[iFreq];
}
double GasInterface::absorptionOpacity_SI(const GasState& gs, size_t iFreq) const
{
	return 100 * (gs._opacityv[iFreq] - gs._scatteringOpacityv[iFreq]);
}

void GasInterface::zeroOpticalProperties(GasState& gs) const
{
	Array zerov(_frequencyv.size());
	gs._emissivityv = zerov;
	gs._opacityv = zerov;
	gs._scatteringOpacityv = zerov;
	gs._temperature = 0;
	gs._ionizedFraction = 0;
}
