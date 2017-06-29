#include "GasInterfaceImpl.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "FreeBound.h"
#include "FreeFree.h"
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
                  _boundBound(make_unique<HydrogenLevels>(make_shared<HydrogenFromFiles>(5),
                                                          frequencyv)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
}

GasInterfaceImpl::GasInterfaceImpl(unique_ptr<NLevel> boundBound, const Array& frequencyv)
                : _frequencyv(frequencyv), _boundBound(move(boundBound))
{
}

GasInterfaceImpl::~GasInterfaceImpl() = default;

void GasInterfaceImpl::solveBalance(GasState& gs, double n, double Tinit,
                                    const Array& specificIntensityv) const
{
#ifndef SILENT
	double isrf = TemplatedUtils::integrate<double>(_frequencyv, specificIntensityv);

	DEBUG("Solving balance under isrf of "
	      << isrf << " erg / s / cm2 / sr = "
	      << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);
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
		   value indicates whether the temperature should increase (net absorption) or
		   decrease (net emission). */
		int counter = 0;
		function<int(double)> evaluateBalance = [&](double logT) -> int {
			counter++;
			s = calculateDensities(n, pow(10., logT), specificIntensityv);
			double netPowerIn = heating(s) - cooling(s);
			DEBUG("Cycle " << counter << ": logT = " << logT
			               << "; netHeating = " << netPowerIn << endl
			               << endl);
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(
		                evaluateBalance, logTinit, 4.e-6, logTmax, logTmin);

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

GasInterfaceImpl::Solution
GasInterfaceImpl::calculateDensities(double n, double T, const Array& specificIntensityv) const
{
	Solution s;
	s.n = n;
	s.T = T;
	s.specificIntensityv = specificIntensityv;

	if (n > 0)
	{
		DEBUG("Calculating state for T = " << T << "K" << endl);

		// Ionization balance
		s.f = Ionization::ionizedFraction(n, T, _frequencyv, specificIntensityv);
		DEBUG("Ionized fraction = " << s.f << endl);

		// Level balance
		double nAtm = n * (1 - s.f);
		double np = n * s.f;
		double ne = np;
		s.levelSolution = _boundBound->solveBalance(nAtm, ne, np, T, specificIntensityv);
	}
	else
	{
		s.f = 0;
		s.levelSolution = _boundBound->solveBalance(0, 0, 0, T, specificIntensityv);
	}
	return s;
}

Array GasInterfaceImpl::emissivityv(const Solution& s) const
{
	const Array& lineEmv = _boundBound->emissivityv(s.levelSolution);
	Array contEmCoeffv(_frequencyv.size());
	_freeBound->addEmissionCoefficientv(s.T, contEmCoeffv);
	_freeFree->addEmissionCoefficientv(s.T, contEmCoeffv);
	return lineEmv + (np_ne(s) / Constant::FPI) * contEmCoeffv;
}

Array GasInterfaceImpl::opacityv(const Solution& s) const
{
	const Array& lineOp = _boundBound->opacityv(s.levelSolution);

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
	return _boundBound->cooling(s.levelSolution);
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	return _boundBound->heating(s.levelSolution);
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

void GasInterfaceImpl::testHeatingCurve(double n, const Array& specificIntensityv) const
{
	const string tab = "\t";
	const int samples = 200;

	double T = 10;
	double factor = pow(1000000. / 10., 1. / samples);

	ofstream output = IOTools::ofstreamFile("heatingcurve.dat");
	output << "# 0temperature 1netHeating 2absorption 3emission"
	       << " 4lineHeating 5lineAbsorption 6lineEmission"
	       << " 7continuumHeating 8continuumAbsorption 9continuumEmission"
	       << " 10ionizationFraction" << endl;
	for (int N = 0; N < samples; N++, T *= factor)
	{
		Solution s = calculateDensities(n, T, specificIntensityv);
		double heat = heating(s);
		double cool = cooling(s);
		double lHeat = lineHeating(s);
		double lCool = lineCooling(s);
		double cHeat = continuumHeating(s);
		double cCool = continuumCooling(s);

		double netHeating = heat - cool;
		double netLine = lHeat - lCool;
		double netCont = cHeat - cCool;

		output << T << tab << netHeating << tab << heat << tab << cool << tab << netLine
		       << tab << lHeat << tab << lCool << tab << netCont << tab << cHeat << tab
		       << cCool << tab << s.f << endl;
	}
	output.close();

	double isrf = TemplatedUtils::integrate<double>(_frequencyv, specificIntensityv);

	DEBUG("Calculated heating curve under isrf of"
	      << isrf << " erg / s / cm2 / sr = "
	      << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);
}
