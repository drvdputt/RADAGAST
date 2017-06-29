#include "GasInterface.h"
#include "Array.h"
#include "GasInterfaceImpl.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "TwoLevelHardcoded.h"

using namespace std;

GasInterface::GasInterface(const valarray<double>& frequencyv)
                : _frequencyv(frequencyv), _pimpl(make_unique<GasInterfaceImpl>(frequencyv))
{
}

GasInterface::GasInterface(const valarray<double>& frequencyv, const std::string& setupChoice)
                : _frequencyv(frequencyv)
{
	unique_ptr<NLevel> boundBound;
	if (!setupChoice.compare("twolevel"))
		boundBound = make_unique<NLevel>(make_shared<TwoLevelHardcoded>(), frequencyv);
	else if (!setupChoice.compare("hhc"))
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenHardcoded>(), frequencyv);
	else if (!setupChoice.compare("hff2"))
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(2), frequencyv);
	else
		boundBound = make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(), frequencyv);
	_pimpl = make_unique<GasInterfaceImpl>(move(boundBound), frequencyv);
}

GasInterface::~GasInterface() = default;

//void GasInterface::setFrequencyv(const Array& frequencyv)
//{
//	_frequencyv = frequencyv;
//	_pimpl = make_unique<GasInterfaceImpl>(frequencyv);
//}

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

void GasInterface::testHeatingCurve(double n, const std::valarray<double>& specificIntensityv) const
{
	_pimpl->testHeatingCurve(n, specificIntensityv);
}
