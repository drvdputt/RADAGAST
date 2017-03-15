#include "HydrogenCalculator.h"
#include "Constants.h"
#include "flags.h"
#include "FreeBound.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "ReadData.h"
#include "TemplatedUtils.h"
#include "TwoLevel.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

HydrogenCalculator::HydrogenCalculator(const vector<double>& frequencyv) :
		_frequencyv(frequencyv), _n(0), _p_specificIntensityv(nullptr), _T(0), _ionizedFraction(0), _levels(
				std::make_unique<TwoLevel>(frequencyv)), _freeBound(std::make_unique<FreeBound>(frequencyv))
{
}

HydrogenCalculator::~HydrogenCalculator() = default;

void HydrogenCalculator::solveBalance(double n, double Tinit, const vector<double>& specificIntensity)
{
	_n = n;
	_p_specificIntensityv = &specificIntensity;

#ifndef SILENT
	double isrf = NumUtils::integrate<double>(_frequencyv, specificIntensity);
#endif
	DEBUG("Solving balance under isrf of " << isrf << " erg / s / cm2 / sr = " << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);

	if (_n > 0)
	{
		const double Tmax = 1000000.;
		const double Tmin = 10;
		const double logTmax = log10(Tmax);
		const double logTmin = log10(Tmin);

		double logTinit = log10(Tinit);

		// Lambda function that will be used by the search algorithm.
		// The state of the system will be updated every time the algorithm calls this function.
		// The return value indicates whether the temperature should increase (net absorption)
		// or decrease (net emission).
		int counter = 0;
		function<int(double)> evaluateBalance = [this, &counter] (double logT) -> int
		{
			counter++;
			calculateDensities(pow(10., logT));
			double netPowerIn = absorption() - emission();
#ifdef VERBOSE
			DEBUG("Cycle " << counter << ": logT = " << logT << "; netHeating = " << netPowerIn << endl << endl);
#endif
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(evaluateBalance, logTinit, 4.e-3,
				logTmax, logTmin);

		// Evaluate the densities for one last time, using the final temperature
		calculateDensities(pow(10., logTfinal));
	}
	else
	{
		calculateDensities(0);
	}
}

void HydrogenCalculator::solveInitialGuess(double n, double T, const vector<double>& p_specificIntensityv)
{
	_n = n;
	_p_specificIntensityv = &p_specificIntensityv;
	calculateDensities(T);
}

GasState HydrogenCalculator::exportState() const
{
	return GasState(_frequencyv, emissivityv(), opacityv(), scatteringOpacityv(), _T, _ionizedFraction);
}

vector<double> HydrogenCalculator::emissivityv() const
{
	vector<double> result;
	const vector<double>& lineEmv = _levels->totalEmissivityv();
	const vector<double>& contCoeffv = _freeBound->emissionCoefficientv(_T);
	for (size_t i = 0; i < _frequencyv.size(); i++)
	{
		double np_ne = _n * _n * _ionizedFraction * _ionizedFraction;
		result.push_back(lineEmv[i] + np_ne * contCoeffv[i] / Constant::FPI);
	}
	return result;
}

vector<double> HydrogenCalculator::opacityv() const
{
	vector<double> result;
	result.reserve(_frequencyv.size());
	const vector<double>& lineOp = _levels->totalOpacityv();
	for (size_t i = 0; i < _frequencyv.size(); i++)
	{
		double ionizOp_i = _n * (1. - _ionizedFraction) * Ionization::crossSection(_frequencyv[i]);
		result.push_back(ionizOp_i + lineOp[i]);
	}
	return result;
}

vector<double> HydrogenCalculator::scatteringOpacityv() const
{
	return _levels->scatteringOpacityv();
}

vector<double> HydrogenCalculator::scatteredv() const
{
	// only the line can scatter
	const vector<double>& opv = _levels->scatteringOpacityv();
	vector<double> intensityOpacityv;
	intensityOpacityv.reserve(opv.size());
	auto op = opv.begin();
	auto opEnd = opv.end();
	auto I_nu = _p_specificIntensityv->begin();
	while (op != opEnd)
	{
		intensityOpacityv.push_back(*I_nu * *op);
		op++;
		I_nu++;
	}
	return intensityOpacityv;
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
	const vector<double>& lineEmv = _levels->totalEmissivityv();
	return Constant::FPI * NumUtils::integrate<double>(_frequencyv, lineEmv);
}

double HydrogenCalculator::lineAbsorption() const
{
	const vector<double>& opv = _levels->totalOpacityv();
	vector<double> intensityOpacityv;
	intensityOpacityv.reserve(opv.size());
	auto op = opv.begin();
	auto opEnd = opv.end();
	auto I_nu = _p_specificIntensityv->begin();
	while (op != opEnd)
	{
		intensityOpacityv.push_back(*I_nu * *op);
		op++;
		I_nu++;
	}
	return Constant::FPI * NumUtils::integrate<double>(_frequencyv, intensityOpacityv);
}

double HydrogenCalculator::continuumEmission() const
{
	const vector<double>& gamma_nu = _freeBound->emissionCoefficientv(_T);
	// emissivity = ne np / 4pi * gamma
	// total emission = 4pi integral(emissivity) = ne np integral(gamma)
	return _n * _n * _ionizedFraction * _ionizedFraction
			* NumUtils::integrate<double>(_frequencyv, gamma_nu);
}

double HydrogenCalculator::continuumAbsorption() const
{
	double atomDensity = _n * (1. - _ionizedFraction);

	vector<double> intensityOpacityv;
	intensityOpacityv.reserve(_frequencyv.size());
	auto I_nu = _p_specificIntensityv->begin();
	for (auto freq = _frequencyv.begin(); freq != _frequencyv.end(); freq++, I_nu++)
		intensityOpacityv.push_back(*I_nu * atomDensity * Ionization::crossSection(*freq));

	return Constant::FPI * NumUtils::integrate<double>(_frequencyv, intensityOpacityv);
}

void HydrogenCalculator::testHeatingCurve()
{
	const string tab = "\t";
	const int samples = 200;

	double T = 10;
	double factor = pow(500000. / 10., 1. / samples);

	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/heating.dat");
	out << "# frequency netHeating absorption emission" << " lineHeating lineAbsorption lineEmission"
			<< " continuumHeating continuumAbsorption continuumEmission" << " ionizationFraction"
			<< endl;
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

		out << T << tab << netHeating << tab << abs << tab << em << tab << netLine << tab << lineAbs
				<< tab << lineEm << tab << netCont << tab << contAbs << tab << contEm << tab
				<< _ionizedFraction << endl;
	}
	out.close();
}

void HydrogenCalculator::calculateDensities(double T)
{
	_T = T;
	if (_n > 0)
	{
#ifdef VERBOSE
		DEBUG("Calculating state for T = " << T << "K" << endl);
#endif

		_ionizedFraction = Ionization::ionizedFraction(_n, T, _frequencyv, *_p_specificIntensityv);
#ifdef VERBOSE
		DEBUG("Ionized fraction = " << _ionizedFraction << endl);
#endif

		double nAtomic = _n * (1 - _ionizedFraction);

		double np = _n * _ionizedFraction;
		double ne = np;
		double alphaTotal = Ionization::recombinationRate(T);

		// approximations from Draine's book, p 138, valid for 3000 to 30000 K
		double T4 = T / 1.e4;
		double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4)); // yes, this is natural log
		double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
		//DEBUG("alphaGround " << alphaGround << " alpha2p " << alpha2p << endl);

		vector<double> sourcev =
		{ ne * np * alphaGround, ne * np * alpha2p };

		// The ionization rate calculation makes no distinction between the levels.
		// When the upper level population is small, and its decay rate is large,
		// the second term doesn't really matter.
		// Therefore, we choose the sink to be the same for each level. Moreover,
		// total source = total sink
		// so we want sink*n0 + sink*n1 = source => sink = totalsource / n because n0/n + n1/n = 1
		double sink = ne * np * alphaTotal / nAtomic;
		vector<double> sinkv =
		{ sink, sink };

		_levels->solveBalance(nAtomic, ne, np, T, *_p_specificIntensityv, sourcev, sinkv);
	}
	else
	{
		_ionizedFraction = 0;
		_levels->solveBalance(0, 0, 0, 0, *_p_specificIntensityv, {0, 0}, {0, 0});
	}
}

