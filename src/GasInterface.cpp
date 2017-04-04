#include "GasInterface.h"
#include "Array.h"
#include "HydrogenCalculator.h"

using namespace std;

GasInterface::GasInterface(const valarray<double>& frequencyv)
                : _frequencyv(frequencyv), _hc(make_unique<HydrogenCalculator>(frequencyv))
{
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gs, double density, double Tinit,
                                  const valarray<double>& specificIntensityv)
{
	gs._previousISRFv = specificIntensityv;
	if (density > 0)
		_hc->solveBalance(gs, density, Tinit, specificIntensityv);
	else
		zeroOpticalProperties(gs);
}

void GasInterface::initializeGasState(GasState& gs, double density, double temperature)
{
	gs._previousISRFv = Array(_frequencyv.size());
	if (density > 0)
		_hc->solveInitialGuess(GasState&, density, temperature);
	else
		zeroOpticalProperties(gs);
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
