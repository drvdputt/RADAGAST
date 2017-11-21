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
	// Choose a level model
	unique_ptr<NLevel> atomicLevels;

	if (atomChoice == "twolevel")
	{
		// Generic level model with hardcoded 2-level data
		auto twoLevelData{make_shared<TwoLevelHardcoded>()};
		atomicLevels = make_unique<NLevel>(twoLevelData, _frequencyv);
	}
	else
	{
		shared_ptr<HydrogenDataProvider> levelData;

		// Choose a level data provider
		if (atomChoice == "hhc")
			levelData = make_shared<HydrogenHardcoded>();
		else if (atomChoice == "hff2")
			levelData = make_shared<HydrogenFromFiles>(2);
		else if (atomChoice == "hff4")
			levelData = make_shared<HydrogenFromFiles>(4);
		else
			levelData = make_shared<HydrogenFromFiles>();

		// Specialized hydrogen level model with the chosen set of level data
		atomicLevels = make_unique<HydrogenLevels>(levelData, _frequencyv);
	}

	// Give this level model to the main implementation.
	_pimpl = make_unique<GasInterfaceImpl>(move(atomicLevels), moleculeChoice, _frequencyv);
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gasState, double n, double Tinit,
                                  const valarray<double>& specificIntensityv,
                                  const GrainInterface& grainInfo) const
{
	gasState._previousISRFv = specificIntensityv;
	if (n > 0)
		_pimpl->solveBalance(gasState, n, Tinit, specificIntensityv, grainInfo);
	else
		zeroOpticalProperties(gasState);
}

void GasInterface::initializeGasState(GasState& gasState, double n, double T,
                                      const GrainInterface& grainInfo) const
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
	gs._densityv = Array{0, 0, 0, 0};
}
} /* namespace GasModule */
