#include "HydrogenCalculator.h"
#include "Constants.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "IonizationBalance.h"
#include "NLevel.h"
#include "NumUtils.h"
#include "TemplatedUtils.h"
#include "Testing.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "global.h"

using namespace std;

HydrogenCalculator::HydrogenCalculator(const Array& frequencyv)
                : _frequencyv(frequencyv), _levels(make_unique<NLevel>(frequencyv)),
                  _freeBound(make_unique<FreeBound>(frequencyv)),
                  _freeFree(make_unique<FreeFree>(frequencyv))
{
	DEBUG("Constructed HydrogenCalculator" << endl);
}

HydrogenCalculator::HydrogenCalculator(Array& frequencyv, bool suggestGrid)
{
	_levels = make_unique<NLevel>(frequencyv);

	int numLines;
	Array lineFreqv, lineWidthv;
	_levels->lineInfo(numLines, lineFreqv, lineWidthv);

	if (suggestGrid)
	{
		double lineWindowFactor = 5;
		vector<double> gridVector(begin(frequencyv), end(frequencyv));
		double thermalFactor = sqrt(Constant::BOLTZMAN * 500000 / Constant::HMASS_CGS) /
		                       Constant::LIGHT;
		lineWidthv = lineWindowFactor * (lineWidthv + lineFreqv * thermalFactor);
		Testing::refineFrequencyGrid(gridVector, 31, 3.,
		                    vector<double>(begin(lineFreqv), end(lineFreqv)),
		                    vector<double>(begin(lineWidthv), end(lineWidthv)));
		frequencyv = Array(gridVector.data(), gridVector.size());
	}
	_frequencyv = frequencyv;
	_levels = make_unique<NLevel>(frequencyv);
	_freeBound = make_unique<FreeBound>(frequencyv);
	_freeFree = make_unique<FreeFree>(frequencyv);
}

HydrogenCalculator::~HydrogenCalculator() {}

void HydrogenCalculator::solveBalance(GasState& gs, double n, double Tinit,
                                      const Array& specificIntensityv)
{
#ifndef SILENT
	double isrf = TemplatedUtils::integrate<double>(_frequencyv, specificIntensityv);
#endif
	DEBUG("Solving balance under isrf of "
	      << isrf << " erg / s / cm2 / sr = "
	      << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);

	double ionizedFraction;

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
			calculateDensities(n, pow(10., logT), specificIntensityv, ionizedFraction);
			double netPowerIn = absorption() - emission();
#ifdef VERBOSE
			DEBUG("Cycle " << counter << ": logT = " << logT
			               << "; netHeating = " << netPowerIn << endl
			               << endl);
#endif
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(
		                evaluateBalance, logTinit, 4.e-3, logTmax, logTmin);

		// Evaluate the densities for one last time, using the final temperature.
		calculateDensities(n, pow(10., logTfinal), specificIntensityv, ionizedFraction);
	}
	else
	{
		calculateDensities(0, 0, specificIntensityv, ionizedFraction);
	}
}

void HydrogenCalculator::solveInitialGuess(GasState& gs, double n, double T)
{
	Array isrfGuess(_frequencyv.size());
	// i'll implement this somewhere else later
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];
		isrfGuess[iFreq] = 2 * Constant::PLANCK * freq * freq * freq / Constant::LIGHT /
		                   Constant::LIGHT /
		                   expm1(Constant::PLANCK * freq / Constant::BOLTZMAN / T);
	}

	// this should be different too. Maybe work with shared pointers or something. Or dump
	// the idea of keeping this in the state of the HC alltogether.
	solveBalance(gs, n, T, isrfGuess);
}

void HydrogenCalculator::calculateDensities(double n, double T, const Array& specificIntensityv, double& ionizedFraction)
{
	if (n > 0)
	{
#ifdef VERBOSE
		DEBUG("Calculating state for T = " << T << "K" << endl);
#endif
		ionizedFraction = Ionization::ionizedFraction(n, T, _frequencyv,
		                                               specificIntensityv);
#ifdef VERBOSE
		DEBUG("Ionized fraction = " << _ionizedFraction << endl);
#endif

		double np = n * ionizedFraction;
		double ne = np;
		double alphaTotal = Ionization::recombinationRate(T);

		// approximations from Draine's book, p 138, valid for 3000 to 30000 K
		// yes, this is natural log
		double T4 = T / 1.e4;
		double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4));
		double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
		double alpha2s = alpha2p / 3.; // TODO: don't forget to change this!

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
		double nAtomic = n * (1. - ionizedFraction);
		double sink = ne * np * alphaTotal / nAtomic;
		Array sinkv({sink, sink, sink, sink, sink, sink});

		_levels->solveBalance(nAtomic, ne, np, T, specificIntensityv, sourcev, sinkv);
	}
	else
	{
		Array zero(_levels->N());
		ionizedFraction = 0;
		_levels->solveBalance(0, 0, 0, T, specificIntensityv, zero, zero);
	}
}

Array HydrogenCalculator::emissivityv() const
{
	const Array& lineEmv = _levels->emissivityv();
	Array contEmCoeffv(_frequencyv.size());
	_freeBound->addEmissionCoefficientv(_temperature, contEmCoeffv);
	_freeFree->addEmissionCoefficientv(_temperature, contEmCoeffv);
	return lineEmv + (np_ne() / Constant::FPI) * contEmCoeffv;
}

Array HydrogenCalculator::opacityv() const
{
	const Array& lineOp = _levels->opacityv();

	Array contOpCoeffv(_frequencyv.size());
	_freeFree->addOpacityCoefficientv(_temperature, contOpCoeffv);
	double npne = np_ne();
	double n_H0 = nAtomic();

	Array totalOp(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double ionizOp_iFreq = n_H0 * Ionization::crossSection(_frequencyv[iFreq]);
		totalOp[iFreq] = ionizOp_iFreq + npne * contOpCoeffv[iFreq] + lineOp[iFreq];
	}
	return totalOp;
}

Array HydrogenCalculator::scatteringOpacityv() const { return _levels->scatteringOpacityv(); }

Array HydrogenCalculator::scatteredv() const
{
	// only the lines can scatter
	const Array& scaOp = _levels->scatteringOpacityv();
	return *_p_specificIntensityv * scaOp;
}

double HydrogenCalculator::emission() const
{
	double lineEm = lineEmission();
	double contEm = continuumEmission();
#ifdef VERBOSE
	DEBUG("emission: line / cont = " << lineEm << " / " << contEm << endl);
#endif
	return lineEm + contEm;
}

double HydrogenCalculator::absorption() const
{
	double lineAbs = lineAbsorption();
	double contAbs = continuumAbsorption();
#ifdef VERBOSE
	DEBUG("absorption: line / cont = " << lineAbs << " / " << contAbs << endl);
#endif
	return lineAbs + contAbs;
}

double HydrogenCalculator::lineEmission() const
{
	return Constant::FPI *
	       TemplatedUtils::integrate<double>(_frequencyv, _levels->emissivityv());
}

double HydrogenCalculator::lineAbsorption() const
{
	Array intensityOpacityv = *_p_specificIntensityv * _levels->opacityv();
	return Constant::FPI * TemplatedUtils::integrate<double>(_frequencyv, intensityOpacityv);
}

double HydrogenCalculator::continuumEmission() const
{
	Array gamma_nuv(_frequencyv.size());
	_freeBound->addEmissionCoefficientv(_temperature, gamma_nuv);
	_freeFree->addEmissionCoefficientv(_temperature, gamma_nuv);

	// emissivity = ne np / 4pi * gamma
	// total emission = 4pi integral(emissivity) = ne np integral(gamma)
	return np_ne() * TemplatedUtils::integrate<double>(_frequencyv, gamma_nuv);
}

double HydrogenCalculator::continuumAbsorption() const
{
	const Array& I_nu = *_p_specificIntensityv;

	double npne = np_ne();
	double nH0 = nAtomic();
	Array freefreeOpCoefv(_frequencyv.size());
	_freeFree->addOpacityCoefficientv(_temperature, freefreeOpCoefv);

	Array intensityOpacityv(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double totalOpacity = nH0 * Ionization::crossSection(_frequencyv[iFreq]) +
		                      npne * freefreeOpCoefv[iFreq];
		intensityOpacityv[iFreq] = I_nu[iFreq] * totalOpacity;
	}
	return Constant::FPI * TemplatedUtils::integrate<double>(_frequencyv, intensityOpacityv);
}
/*
void HydrogenCalculator::testHeatingCurve()
{
	const string tab = "\t";
	const int samples = 200;

	double T = 10;
	double factor = pow(500000. / 10., 1. / samples);

	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/heating.dat");
	out << "# 0frequency 1netHeating 2absorption 3emission"
	    << " 4lineHeating 5lineAbsorption 6lineEmission"
	    << " 7continuumHeating 8continuumAbsorption 9continuumEmission"
	    << " 10ionizationFraction" << endl;
	for (int s = 0; s < samples; s++, T *= factor)
	{
		calculateDensities(T);
		double abs = absorption();
		double em = emission();
		double lineAbs = lineAbsorption();
		double lineEm = lineEmission();
		double contAbs = continuumAbsorption();
		double contEm = continuumEmission();

		double netHeating = abs - em;
		double netLine = lineAbs - lineEm;
		double netCont = contAbs - contEm;

		out << T << tab << netHeating << tab << abs << tab << em << tab << netLine << tab
		    << lineAbs << tab << lineEm << tab << netCont << tab << contAbs << tab << contEm
		    << tab << _ionizedFraction << endl;
	}
	out.close();
}
*/
