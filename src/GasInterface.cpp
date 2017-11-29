#include "GasInterface.h"
#include "Array.h"
#include "Error.h"
#include "GasInterfaceImpl.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "TwoLevelHardcoded.h"

#include <sstream>

using namespace std;

namespace GasModule
{
GasInterface::GasInterface(const valarray<double>& frequencyv, const string& atomChoice,
                           const string& moleculeChoice)
                : _frequencyv(frequencyv)
{
	// Level model for the atom.
	unique_ptr<NLevel> atomicLevels;
	// Choose a data provider and make a level model out of it
	if (atomChoice == "twolevel")
	{
		// Generic level model with hardcoded 2-level data
		auto twoLevelData{make_shared<TwoLevelHardcoded>()};
		atomicLevels = make_unique<NLevel>(twoLevelData, _frequencyv);
	}
	else
	{
		// Specialized level model for atomic hydrogen
		shared_ptr<HydrogenDataProvider> hLevelData;

		// Choose a hydrogen data provider
		if (atomChoice == "hhc")
			hLevelData = make_shared<HydrogenHardcoded>();
		else if (atomChoice == "hff2")
			hLevelData = make_shared<HydrogenFromFiles>(2);
		else if (atomChoice == "hff4")
			hLevelData = make_shared<HydrogenFromFiles>(4);
		else
			hLevelData = make_shared<HydrogenFromFiles>();

		atomicLevels = make_unique<HydrogenLevels>(hLevelData, _frequencyv);
	}

	// Level model for the molecule.
	unique_ptr<H2Levels> molecularLevels;
	if (moleculeChoice == "none")
		molecularLevels = nullptr;
	else
	{
		// Create a data provider.
		shared_ptr<H2FromFiles> h2LevelData;

		if (moleculeChoice.empty())
			h2LevelData = make_shared<H2FromFiles>();
		else
		{
			int maxJ, maxV;
			istringstream(moleculeChoice) >> maxJ >> maxV;
			if (maxJ < 0 || maxV < 0)
				Error::runtime("moleculeChoice is not of a correct format. It "
				               "should be "
				               "\"maxJ maxV\"");
			h2LevelData = make_shared<H2FromFiles>(maxJ, maxV);
		}
		// Make a level model out of it.
		molecularLevels = make_unique<H2Levels>(h2LevelData, _frequencyv);
	}

	// Give these level models to the main implementation (TODO: check if giving a nullptr works
	// correctly)
	_pimpl = make_unique<GasInterfaceImpl>(move(atomicLevels), move(molecularLevels),
	                                       _frequencyv);
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
