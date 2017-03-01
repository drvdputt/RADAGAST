#include "Constants.h"
#include "GasSpecies.h"
#include "NumUtils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>
#include <algorithm>
#include "IonizationBalance.h"
#include "ReadData.h"

namespace
{

template<typename T>
T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax, T xMin)
{
	if (xInit > xMax || xInit < xMin)
		throw "xInit is out of given range for binary search";

	double upperbound = xMax;
	double lowerbound = xMin;
	double current = xInit;

	while (upperbound - lowerbound > xTolerance)
	{
		// Indicates whether the (unknown) final value lies above or below the current one
		double iCompare = searchDirection(current);

		if (iCompare > 0)
			lowerbound = current;
		else if (iCompare < 0)
			upperbound = current;
		else
			return current;

		current = (lowerbound + upperbound) / 2.;
	}
	return current;
}
}

using namespace std;

GasSpecies::GasSpecies(const vector<double>& frequencyv) :
		_frequencyv(frequencyv), _n(0), _p_specificIntensityv(nullptr), _T(0), _ionizedFraction(0), _levels(
				frequencyv)
{
	// Read the free-bound continuum data from Ercolano and Storey (2006)
	// Adapted from NEBULAR source code by M. Schirmer (2016)
	vector<double> fileFrequencyv;
	vector<vector<double>> fileGammaDaggervv;

	string file("/Users/drvdputt/GasModule/git/dat/t3_elec_reformat.ascii");
	ReadData::recombinationContinuum(file, fileFrequencyv, _thresholdv, _logTemperaturev,
			fileGammaDaggervv);

	size_t numcol = _logTemperaturev.size();
	size_t numrow = fileFrequencyv.size();

	// Now resample this data according to the frequency grid
	_gammaDaggervv.resize(_frequencyv.size(), vector<double>(numcol, 0));

	cout << "frequency range from file: " << fileFrequencyv[0] << " to " << fileFrequencyv.back() << endl;
	cout << "frequency range: " << _frequencyv[0] << " to " << _frequencyv.back() << endl;

	// Then, apply a linear interpolation across the frequencies (rows) for every temperature (column)
	for (size_t col = 0; col < numcol; col++)
	{
		// Extract the column
		vector<double> column;
		column.reserve(numrow);
		for (size_t row = 0; row < numrow; row++)
			column.push_back(fileGammaDaggervv[row][col]);

		// Resample it
		const vector<double>& column_resampled = NumUtils::interpol<double>(column, fileFrequencyv,
				frequencyv, -1, -1);

		// And copy it over
		for (size_t row = 0; row < _gammaDaggervv.size(); row++)
			_gammaDaggervv[row][col] = column_resampled[row];
	}

//#define _PRINT_CONTINUUM_DATA_
#ifdef _PRINT_CONTINUUM_DATA_
	// DEBUG: print out the table as read from the file
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/loadedContinuum.dat");
	for (size_t iNu = 0; iNu < fileGammaDaggervv.size(); iNu++)
	{
		for (double d : fileGammaDaggervv[iNu])
		out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: print out the interpolated table
	out.open("/Users/drvdputt/GasModule/run/interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(); iNu++)
	{
		for (double d : _gammaDaggervv[iNu])
		out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: Test the temperature interpolation function (at least a copy pasta of a part)
	out.open("/Users/drvdputt/GasModule/run/bi-interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(); iNu++)
	{
		for (double logT = 2; logT < 5; logT += 0.01)
		{
			// Find the grid point to the right of the requested log-temperature
			size_t iRight = NumUtils::index(logT, _logTemperaturev);
			// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
			double wRight = (logT - _logTemperaturev[iRight - 1])
			/ (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

			// Interpolate gamma^dagger linearly in log T space
			double gammaDagger = (_gammaDaggervv[iNu][iRight - 1] * (1 - wRight)
					+ _gammaDaggervv[iNu][iRight] * wRight);
			out << scientific << gammaDagger << "\t";
		}
		out << endl;
	}
	out.close();
#endif /*_PRINT_CONTINUUM_DATA_*/
}

void GasSpecies::solveBalance(double n, double Tinit, const vector<double>& specificIntensity)
{
	const double Tmax = 1000000.;
	const double Tmin = 10;
	const double logTmax = log10(Tmax);
	const double logTmin = log10(Tmin);

	_n = n;
	_p_specificIntensityv = &specificIntensity;

	// Do a binary search in logT space
	double logTinit = log10(Tinit);

	// Lambda function that will be used by the search algorithm.
	// The state of the system will be updated every time the algorithm calls this function.
	// The return value indicates whether the temperature should increase (net absorption)
	// or decrease (net emission).
	int counter;
	function<int(double)> evaluateBalance = [this, &counter] (double logT) -> int
	{
		counter++;
		calculateDensities(pow(10., logT));
		double netPowerIn = absorption() - emission();
		cout << "Cycle " << counter << ": logT = " << logT << "; netHeating = "
		<< netPowerIn << endl << endl;
		return (netPowerIn > 0) - (netPowerIn < 0);
	};

	double logTfinal = binaryIntervalSearch<double>(evaluateBalance, logTinit, 4.e-3, logTmax, logTmin);

	// Evaluate the densities for one last time, using the final temperature
	calculateDensities(pow(10., logTfinal));
}

vector<double> GasSpecies::emissivityv() const
{
	vector<double> result;
	const vector<double>& lineEmv = _levels.emissivityv();
	const vector<double>& contCoeffv = recombinationEmissionCoeff(_T);
	for (size_t i = 0; i < _frequencyv.size(); i++)
	{
		double np_ne = _n * _n * _ionizedFraction * _ionizedFraction;
		result.push_back(lineEmv[i] + np_ne * contCoeffv[i] / Constant::FPI);
	}
	return result;
}

vector<double> GasSpecies::opacityv() const
{
	vector<double> result;
	result.reserve(_frequencyv.size());
	const vector<double>& lineOp = _levels.opacityv();
	for (size_t i = 0; i < _frequencyv.size(); i++)
	{
		double ionizOp_i = _n * (1. - _ionizedFraction) * Ionization::crossSection(_frequencyv[i]);
		result.push_back(ionizOp_i + lineOp[i]);
	}
	return result;
}

double GasSpecies::emission() const
{
	double lineEm = lineEmission();
	double contEm = continuumEmission();
	cout << "emission: line / cont = " << lineEm << " / " << contEm << endl;
	return lineEm + contEm;
}

double GasSpecies::absorption() const
{
	double lineAbs = lineAbsorption();
	double contAbs = continuumAbsorption();
	cout << "absorption: line / cont = " << lineAbs << " / " << contAbs << endl;
	return lineAbs + contAbs;
}

double GasSpecies::lineEmission() const
{
	const vector<double>& lineEmv = _levels.emissivityv();
	return Constant::FPI * NumUtils::integrate<double>(_frequencyv, lineEmv);
}

double GasSpecies::lineAbsorption() const
{
	const vector<double>& opv = _levels.opacityv();
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

double GasSpecies::continuumEmission() const
{
	const vector<double>& gamma_nu = recombinationEmissionCoeff(_T);
	// emissivity = ne np / 4pi * gamma
	// total emission = 4pi integral(emissivity) = ne np integral(gamma)
	return _n * _n * _ionizedFraction * _ionizedFraction
			* NumUtils::integrate<double>(_frequencyv, gamma_nu);
}

double GasSpecies::continuumAbsorption() const
{
	double atomDensity = _n * (1. - _ionizedFraction);

	vector<double> intensityOpacityv;
	intensityOpacityv.reserve(_frequencyv.size());
	auto I_nu = _p_specificIntensityv->begin();
	for (auto freq = _frequencyv.begin(); freq != _frequencyv.end(); freq++, I_nu++)
		intensityOpacityv.push_back(*I_nu * atomDensity * Ionization::crossSection(*freq));

	return Constant::FPI * NumUtils::integrate<double>(_frequencyv, intensityOpacityv);
}

void GasSpecies::testHeatingCurve()
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

void GasSpecies::calculateDensities(double T)
{
	_T = T;
	cout << "Calculating state for T = " << T << "K" << endl;

	_ionizedFraction = Ionization::ionizedFraction(_n, T, _frequencyv, *_p_specificIntensityv);

	cout << "Ionized fraction = " << _ionizedFraction << endl;

	double nAtomic = _n * (1 - _ionizedFraction);
	// For now we can sum over the collision partners, as TwoLevel threats all collision partners in the same way
	// neutral + ion + electron densities
	double nTotal = (1 + _ionizedFraction) * _n;

	double ne_np = _n * _n * _ionizedFraction * _ionizedFraction;
	double alphaTotal = Ionization::recombinationRate(T);

	// approximations from Draine's book, p 138, valid for 3000 to 30000 K
	double T4 = T / 1.e4;
	double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4)); // yes, this is natural log
	double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
	//cout << "alphaGround " << alphaGround << " alpha2p " << alpha2p << endl;

	vector<double> sourcev =
	{ ne_np * alphaGround, ne_np * alpha2p };

	// The ionization rate calculation makes no distinction between the levels.
	// When the upper level population is small, and its decay rate is large, the second term doesn't really matter.
	// Therefore, we choose the sink to be the same for each level. Moreover, total source = total sink
	// so we want sink*n0 + sink*n1 = source => sink = totalsource / n because n0/n + n1/n = 1
	double sink = ne_np * alphaTotal / nAtomic;
	vector<double> sinkv =
	{ sink, sink };

	_levels.doLevels(nAtomic, nTotal, T, *_p_specificIntensityv, sourcev, sinkv);
}

vector<double> GasSpecies::recombinationEmissionCoeff(double T) const
{
	double logT = log10(T);

	if (logT > _logTemperaturev.back() || logT < _logTemperaturev[0])
	{
		cout << "Warning: temperature " << T << "K is outside of data range for free-bound continuum"
				<< endl;
		return vector<double>(_frequencyv.size(), 0.);
	}

	vector<double> result;
	result.reserve(_frequencyv.size());

	// Find the grid point to the right of the requested log-temperature
	size_t iRight = NumUtils::index(logT, _logTemperaturev);
	// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
	double wRight = (logT - _logTemperaturev[iRight - 1])
			/ (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

	// We will use equation (1) of Ercolano and Storey 2006 to remove the normalization of the data
	double Ttothe3_2 = pow(T, 3. / 2.);
	double kT = Constant::BOLTZMAN * T;
	// "the nearest threshold of lower energy"
	// i.e. the frequency in _thresholdv lower than the current one
	size_t iThreshold = 0;
	double tE = 0;

	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/gammanu.dat");

	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];

		// Interpolate gamma^dagger linearly in log T space
		double gammaDagger = (_gammaDaggervv[iFreq][iRight - 1] * (1 - wRight)
				+ _gammaDaggervv[iFreq][iRight] * wRight);

		// Skip over zero data, or when we are below the first threshold
		if (!gammaDagger || freq < _thresholdv[0])
		{
			result.push_back(0.);
		}
		else
		{
			// If the last threshold hasn't been passed yet, check if we have passed the next (i + 1)
			if (freq < _thresholdv.back())
			{
				// This block must only be executed when a new threshold is passed
				if (freq > _thresholdv[iThreshold + 1])
				{
					// find the next threshold of lower frequency
					// (don't just pick the next one, as this wouldn't work with very coarse grids)
					iThreshold = NumUtils::index<double>(freq, _thresholdv) - 1;
					tE = Constant::PLANCK * _thresholdv[iThreshold];
				}
			}
			// When we have just moved past the last threshold, set the index one last time
			else if (iThreshold < _thresholdv.size() - 1)
			{
				iThreshold = _thresholdv.size() - 1;
				tE = Constant::PLANCK * _thresholdv.back();
			}
			double E = Constant::PLANCK * freq;

			double normalizationFactor = 1.e34 * Ttothe3_2 * exp((E - tE) / kT);
			double gammaNu = gammaDagger / normalizationFactor;
			out << freq << "\t" << gammaNu / 1.e-40 << endl;
			result.push_back(gammaNu);
		}
	}
	out.close();
	return result;
}

