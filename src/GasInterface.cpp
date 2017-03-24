#include "GasInterface.h"
#include "Array.h"
#include "HydrogenCalculator.h"

using namespace std;

GasInterface::GasInterface(const valarray<double>& frequencyv) :
		_frequencyv(frequencyv), _hc(make_unique < HydrogenCalculator > (frequencyv))
{
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gs, double density, double Tinit,
		const valarray<double>& specificIntensityv)
{
	gs._previousISRFv = specificIntensityv;
	if (density > 0)
	{
		_hc->solveBalance(gs, density, Tinit, specificIntensityv);
		fillOpticalProperties(gs);
	}
	else
	{
		zeroOpticalProperties(gs);
	}
}

void GasInterface::initializeGasState(GasState& gs, double density, double temperature)
{
	gs._previousISRFv = Array(_frequencyv.size());
	if (density > 0)
	{
		_hc->solveInitialGuess(density, temperature);
		fillOpticalProperties(gs);
	}
	else
	{
		zeroOpticalProperties(gs);
	}
}

void GasInterface::fillOpticalProperties(GasState& gs) const
{
	gs._emissivityv = _hc->emissivityv();
	gs._opacityv = _hc->opacityv();
	gs._scatteringOpacityv = _hc->scatteringOpacityv();
	gs._temperature = _hc->temperature();
	gs._ionizedFraction = _hc->ionizedFraction();
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
