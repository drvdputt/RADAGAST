#include "GasInterfaceImpl.h"
#include "Chemistry.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "Testing.h"

using namespace std;

GasInterfaceImpl::GasInterfaceImpl(const Array& frequencyv)
                : _frequencyv(frequencyv),
                  _atomicLevels(make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(5),
                                                            frequencyv)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
}

GasInterfaceImpl::GasInterfaceImpl(unique_ptr<NLevel> atomModel, const Array& frequencyv)
                : _frequencyv(frequencyv), _atomicLevels(std::move(atomModel)),
                  _molecularLevels(make_unique<H2Levels>(make_shared<H2FromFiles>(), frequencyv)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
}

GasInterfaceImpl::~GasInterfaceImpl() = default;

void GasInterfaceImpl::solveInitialGuess(GasState& gs, double n, double T) const
{
	Array isrfGuess(_frequencyv.size());
	// i'll put this somewhere else later
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];
		isrfGuess[iFreq] = SpecialFunctions::planck(freq, T);
	}
	solveBalance(gs, n, T, isrfGuess);
}

void GasInterfaceImpl::solveBalance(GasState& gs, double n, double Tinit,
                                    const Array& specificIntensityv) const
{
#ifndef SILENT
	double isrf = TemplatedUtils::integrate<double>(_frequencyv, specificIntensityv);

	DEBUG("Solving balance under isrf of " << isrf << " erg / s / cm2 / sr = ");
	DEBUG(isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);
#endif
	Solution s;

	if (n > 0)
	{
		const double Tmax = 1000000.;
		const double Tmin = 10;
		const double logTmax = log10(Tmax);
		const double logTmin = log10(Tmin);

		double logTinit = log10(Tinit);

		/* Lambda function that will be used by the search algorithm. The state of the
		   system will be updated every time the algorithm calls this function. The return
		   value indicates whether the temperature should increase (there is net heating so
		   we need a higher temperature leading to more cooling) or decrease (there is net
		   cooling so we need a lower temperature leading to less cooling). */
		int counter = 0;
		function<int(double)> evaluateThermalBalance = [&](double logT) -> int {
			counter++;
			s = calculateDensities(n, pow(10., logT), specificIntensityv);
			double netPowerIn = heating(s) - cooling(s);
			DEBUG("Cycle " << counter << ": logT = " << logT
			               << "; netHeating = " << netPowerIn << endl
			               << endl);
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(
		                evaluateThermalBalance, logTinit, 4.e-6, logTmax, logTmin);

		// Evaluate the densities for one last time, using the final temperature.
		s = calculateDensities(n, pow(10., logTfinal), specificIntensityv);
	}
	else
	{
		s = calculateDensities(0, 0, specificIntensityv);
	}

	const Array& emv = emissivityv(s);
	const Array& opv = opacityv(s);
	const Array& scv = scatteringOpacityv(s);

#ifdef SANITY
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		if (specificIntensityv[iFreq] < 0 || emv[iFreq] < 0 || opv[iFreq] < 0 ||
		    scv[iFreq] < 0)
		{
			cout << "GasModule: negative value in one of the optical properties";
		}
	}
#endif
	// Put the relevant data into the gas state
	gs = GasState(specificIntensityv, emv, opv, scv, s.T, s.f);
}

GasInterfaceImpl::Solution
GasInterfaceImpl::calculateDensities(double n, double T, const Array& specificIntensityv) const
{
	Solution s;
	s.n = n;
	s.T = T;
	s.specificIntensityv = specificIntensityv;

	// Initial guess for the chemistry
	s.chemistrySolution = Eigen::VectorXd(4);
	s.chemistrySolution << 0, 0, n, 0;

	if (n > 0)
	{
		DEBUG("Calculating state for T = " << T << "K" << endl);

		// TODO: Chemical network goes here. Should at least include ne np nH nH2, and will
		// make use of a variety of functions from Ionization, NLevel and others to calculate the
		// reaction rates.
		bool stopCriterion = false;
		while (!stopCriterion)
		{
			s.chemistrySolution = _chemistry->solveBalance(
					reactionRates(T, specificIntensityv, s.H2Solution),
					s.chemistrySolution);

			// TODO: keep this for a while. I will later compare the result of the chemical network to
			// this
			// Ionization balance
			s.f = Ionization::solveBalance(n, T, _frequencyv, specificIntensityv);
			DEBUG("Ionized fraction = " << s.f << endl);

			// Level balance
			double nAtm = n * (1 - s.f);
			double np = n * s.f;
			double ne = np;
			s.HSolution = _atomicLevels->solveBalance(s.chemistrySolution[inH], ne, np, T,
								  specificIntensityv);
			s.H2Solution = _molecularLevels->solveBalance(s.chemistrySolution[inH2], ne, np, T,
								      specificIntensityv);

			// TODO: stopCriterion = some evaluation
		}
	}
	else
	{
		s.f = 0;
		s.HSolution = _atomicLevels->solveBalance(0, 0, 0, T, specificIntensityv);
	}
	return s;
}

Array GasInterfaceImpl::emissivityv(const Solution& s) const
{
	const Array& lineEmv = _atomicLevels->emissivityv(s.HSolution);
	Array contEmCoeffv(_frequencyv.size());
	_freeBound->addEmissionCoefficientv(s.T, contEmCoeffv);
	_freeFree->addEmissionCoefficientv(s.T, contEmCoeffv);
	return lineEmv + (np_ne(s) / Constant::FPI) * contEmCoeffv;
}

Array GasInterfaceImpl::opacityv(const Solution& s) const
{
	const Array& lineOp = _atomicLevels->opacityv(s.HSolution);

	Array contOpCoeffv(_frequencyv.size());
	_freeFree->addOpacityCoefficientv(s.T, contOpCoeffv);
	double npne = np_ne(s);
	double n_H0 = nAtomic(s);

	Array totalOp(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double ionizOp_iFreq = n_H0 * Ionization::crossSection(_frequencyv[iFreq]);
		totalOp[iFreq] = ionizOp_iFreq + npne * contOpCoeffv[iFreq] + lineOp[iFreq];
#ifdef SANITY
		if (totalOp[iFreq] < 0)
		{
			cout << "Negative opacity!";
		}
#endif
	}
	return totalOp;
}

Array GasInterfaceImpl::scatteringOpacityv(const Solution& s) const
{
	return Array(s.specificIntensityv.size());
}

double GasInterfaceImpl::cooling(const Solution& s) const
{
	double lineCool = lineCooling(s);
	double contCool = continuumCooling(s);
	DEBUG("cooling: line / cont = " << lineCool << " / " << contCool << endl);
	return lineCool + contCool;
}

double GasInterfaceImpl::heating(const Solution& s) const
{
	double lineHeat = lineHeating(s);
	double contHeat = continuumHeating(s);
	DEBUG("heating: line / cont = " << lineHeat << " / " << contHeat << endl);
	return lineHeat + contHeat;
}

double GasInterfaceImpl::lineCooling(const Solution& s) const
{
	return _atomicLevels->cooling(s.HSolution);
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	return _atomicLevels->heating(s.HSolution);
}

double GasInterfaceImpl::continuumCooling(const Solution& s) const
{
	return _freeFree->cooling(np_ne(s), s.T) + Ionization::cooling(s.n, s.f, s.T);
}

double GasInterfaceImpl::continuumHeating(const Solution& s) const
{
	return _freeFree->heating(np_ne(s), s.T, s.specificIntensityv) +
	       Ionization::heating(s.n, s.f, s.T, _frequencyv, s.specificIntensityv);
}
