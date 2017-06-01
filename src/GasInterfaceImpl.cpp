#include "GasInterfaceImpl.h"

#include "Constants.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "HydrogenLevels.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "TemplatedUtils.h"
#include "Testing.h"
#include "global.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

GasInterfaceImpl::GasInterfaceImpl(const Array& frequencyv)
                : _frequencyv(frequencyv), _boundBound(make_unique<HydrogenLevels>(frequencyv)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
	DEBUG("Constructed HydrogenCalculator" << endl);
}

/* there's a lot of double work happening here, but this constructor shouldn't be called too often
 */
GasInterfaceImpl::GasInterfaceImpl(const Array& frequencyv, bool improveGrid)
                : _boundBound(make_unique<HydrogenLevels>()),
                  _freeBound(make_unique<FreeBound>(frequencyv))
{
	if (improveGrid)
	{
		// Add extra points for the lines
		int numLines;
		Array lineFreqv, lineWidthv;
		_boundBound->lineInfo(numLines, lineFreqv, lineWidthv);

		double lineWindowFactor = 1.;
		double thermalFactor = sqrt(Constant::BOLTZMAN * 500000 / Constant::HMASS_CGS) /
		                       Constant::LIGHT;
		lineWidthv = lineWindowFactor * (lineWidthv + lineFreqv * thermalFactor);

		vector<double> gridVector(begin(frequencyv), end(frequencyv));
		Testing::refineFrequencyGrid(gridVector, 13, 2.5, lineFreqv, lineWidthv);

		// And for the jumps in the bound-bound spectrum
		const Array& thresholdv = _freeBound->thresholdv();
		// Don't bother with the last jump, because that's the end of the data
		Array jumpFreqv(&thresholdv[0], thresholdv.size() - 2);
		Testing::refineFrequencyGrid(gridVector, 3, 1., jumpFreqv, 1e-6 * jumpFreqv);

		// And for the ionization threshold
		Array ionThr({Ionization::THRESHOLD});
		Testing::refineFrequencyGrid(gridVector, 3, 1., ionThr, 1e-6 * ionThr);

		// Overwrite the frequencygrid with the improved one
		_frequencyv = Array(gridVector.data(), gridVector.size());
	}
	else
		_frequencyv = frequencyv;

	_boundBound->setFrequencyv(_frequencyv);

	// Overwrite these objects with new ones
	_freeBound = make_unique<FreeBound>(_frequencyv);
	_freeFree = make_unique<FreeFree>(_frequencyv);
}

GasInterfaceImpl::~GasInterfaceImpl() {}

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
	gs = GasState(_frequencyv, specificIntensityv, emv, opv, scv, s.T, s.f);
}

void GasInterfaceImpl::solveInitialGuess(GasState& gs, double n, double T) const
{
	Array isrfGuess(_frequencyv.size());
	// i'll put this somewhere else later
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];
		isrfGuess[iFreq] = 2 * Constant::PLANCK * freq * freq * freq / Constant::LIGHT /
		                   Constant::LIGHT /
		                   expm1(Constant::PLANCK * freq / Constant::BOLTZMAN / T);
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
		s.f = Ionization::ionizedFraction(n, T, _frequencyv, specificIntensityv);
		DEBUG("Ionized fraction = " << s.f << endl);

		double np = n * s.f;
		double ne = np;
		double alphaTotal = Ionization::recombinationRateCoeff(T);

		// approximations from Draine's book, p 138, valid for 3000 to 30000 K
		// yes, this is natural log
		double T4 = T / 1.e4;
		double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4));
		double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
		double alpha2s = 2.34e-14 * pow(T4, -0.537 - 0.019 * log(T4));

		// 2015-Raga (A13)
		double t = log10(T4);
		vector<double> logAlpha3poly = {-13.3377, -0.7161, -0.1435, -0.0386, 0.0077};
		vector<double> logAlpha4poly = {-13.5225, -0.7928, -0.1749, -0.0412, 0.0154};
		vector<double> logAlpha5poly = {-13.6820, -0.8629, -0.1957, -0.0375, 0.0199};

		double alpha3 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha3poly));
		double alpha4 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha4poly));
		double alpha5 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha5poly));

		Array sourcev({alphaGround, alpha2p, alpha2s, alpha3, alpha4, alpha5});
		sourcev *= ne * np;

		/* The ionization rate calculation makes no distinction between the levels. When the
		 upper level population is small, and its decay rate is large, the second term
		 doesn't really matter. Therefore, we choose the sink to be the same for each level.
		 Moreover, total source = total sink so we want sink*n0 + sink*n1 = source => sink =
		 totalsource / n because n0/n + n1/n = 1.
		 */
		double nAtm = n * (1. - s.f);
		double sink = ne * np * alphaTotal / nAtm;
		Array sinkv({sink, sink, sink, sink, sink, sink});

		s.levelSolution = _boundBound->solveBalance(nAtm, ne, np, T, specificIntensityv,
		                                            sourcev, sinkv);
	}
	else
	{
		s.f = 0;
		Array zero(_boundBound->nLv());
		s.levelSolution = _boundBound->solveBalance(0, 0, 0, T, specificIntensityv, zero,
		                                            zero);
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
	//	return Constant::FPI *
	//	       TemplatedUtils::integrate<double>(_frequencyv,
	//	                                         _boundBound->emissivityv(s.levelSolution));
	return _boundBound->cooling(s.levelSolution);
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	//	Array intensityOpacityv = s.specificIntensityv *
	//_boundBound->opacityv(s.levelSolution); 	return Constant::FPI *
	// TemplatedUtils::integrate<double>(_frequencyv, intensityOpacityv);
	return _boundBound->heating(s.levelSolution);
}

double GasInterfaceImpl::continuumCooling(const Solution& s) const
{
	+	       Ionization::cooling(s.n, s.f, s.T);
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

	ofstream output;
	output.open("/Users/drvdputt/GasModule/run/heating.dat");
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
