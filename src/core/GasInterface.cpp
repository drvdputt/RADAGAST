#include "GasInterface.h"
#include "Array.h"
#include "Error.h"
#include "GasInterfaceImpl.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "Timer.h"
#include "TwoLevelHardcoded.h"

#include <sstream>

using namespace std;

namespace GasModule
{
GasInterface::GasInterface(const valarray<double>& iFrequencyv,
                           const valarray<double>& oFrequencyv,
                           const valarray<double>& eFrequencyv,
                           const string& atomChoice, const string& moleculeChoice)
                : _iFrequencyv{iFrequencyv}, _oFrequencyv{oFrequencyv},
                  _eFrequencyv{eFrequencyv}
{
	unique_ptr<HydrogenLevels> atomicLevels;

	// Choose a hydrogen data provider
	shared_ptr<HydrogenDataProvider> hLevelData;
	if (atomChoice == "hhc")
		hLevelData = make_shared<HydrogenHardcoded>();
	else if (atomChoice == "hff2")
		hLevelData = make_shared<HydrogenFromFiles>(2);
	else if (atomChoice == "hff4")
		hLevelData = make_shared<HydrogenFromFiles>(4);
	else
		hLevelData = make_shared<HydrogenFromFiles>();

	atomicLevels = make_unique<HydrogenLevels>(hLevelData);

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
		molecularLevels = make_unique<H2Levels>(h2LevelData);
	}

	// Give these level models to the main implementation
	_pimpl = make_unique<GasInterfaceImpl>(move(atomicLevels), move(molecularLevels));
}

GasInterface::~GasInterface() = default;

void GasInterface::updateGasState(GasState& gasState, double n, double Tinit,
                                  const valarray<double>& specificIntensityv,
                                  const GrainInterface& grainInfo, GasDiagnostics* gd) const
{
	Timer t("Update gas state");
	// Create a spectrum object which makes it easier to pass around the frequencies and the
	// values for the specific intensity. It can also be used to easily interpolate the
	// specific intensity for a certain frequency.
	Spectrum specificIntensity(_iFrequencyv, specificIntensityv);
	if (n > 0)
	{
		GasInterfaceImpl::Solution s = _pimpl->solveTemperature(n, Tinit, specificIntensity, grainInfo);
		// To fix the temperature, use this:
		// GasInterfaceImpl::Solution s = _pimpl->solveDensities(n, 6000., specificIntensity, grainInfo);
		gasState = _pimpl->makeGasStateFromSolution(s, _oFrequencyv, _eFrequencyv);
		if (gd)
			_pimpl->fillGasDiagnosticsFromSolution(s, grainInfo, gd);
	}
	else
		zeroOpticalProperties(gasState);
}

void GasInterface::initializeGasState(GasState& gasState, double n, double T,
                                      const GrainInterface& grainInfo, GasDiagnostics* gd) const
{
	if (n > 0)
	{
		GasInterfaceImpl::Solution s = _pimpl->solveInitialGuess(n, T, grainInfo);
		gasState = _pimpl->makeGasStateFromSolution(s, _oFrequencyv, _eFrequencyv);
		if (gd)
			_pimpl->fillGasDiagnosticsFromSolution(s, grainInfo, gd);
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
