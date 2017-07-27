#include "GasInterface.h"
#include "Array.h"
#include "GasInterfaceImpl.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "TwoLevelHardcoded.h"

using namespace std;

namespace GasModule
{
GasInterface::GasInterface(const valarray<double>& frequencyv, const string& atomChoice,
                           bool moleculeChoice)
                : _frequencyv(frequencyv)
{
	unique_ptr<NLevel> boundBound;
	if (!atomChoice.compare("twolevel"))
		boundBound = make_unique<NLevel>(make_shared<TwoLevelHardcoded>(), frequencyv);
	else if (!atomChoice.compare("hhc"))
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenHardcoded>(),
		                                         frequencyv);
	else if (!atomChoice.compare("hff2"))
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(2),
		                                         frequencyv);
	else if (!atomChoice.compare("hff4"))
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(4),
		                                         frequencyv);
	else
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(),
		                                         frequencyv);
	_pimpl = make_unique<GasInterfaceImpl>(move(boundBound), moleculeChoice, frequencyv);
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gasState, double n, double Tinit,
                                  const valarray<double>& specificIntensityv,
                                  const GrainInfo& grainInfo) const
{
	gasState._previousISRFv = specificIntensityv;
	if (n > 0)
		_pimpl->solveBalance(gasState, n, Tinit, specificIntensityv, grainInfo);
	else
		zeroOpticalProperties(gasState);
}

void GasInterface::initializeGasState(GasState& gasState, double n, double T,
                                      const GrainInfo& grainInfo) const
{
	gasState._previousISRFv = Array(_frequencyv.size());
	if (n > 0)
		_pimpl->solveInitialGuess(gasState, n, T, grainInfo);
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
	Array zerov(_frequencyv.size());
	gs._emissivityv = zerov;
	gs._opacityv = zerov;
	gs._scatteringOpacityv = zerov;
	gs._temperature = 0;
	gs._ionizedFraction = 0;
}
} /* namespace GasModule */
