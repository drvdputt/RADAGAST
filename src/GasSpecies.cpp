#include "Constants.h"
#include "GasSpecies.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>
#include <algorithm>

using namespace std;

GasSpecies::GasSpecies(const vector<double>& wavelengthv)
	: _wavelengthv(wavelengthv), _n(0), _T(0), _ionizedFraction(0), _levels(wavelengthv)
{
	// Read the free-bound continuum data from Ercolano and Storey (2006)
	// Adapted from NEBULAR source code by M. Schirmer (2016)
	size_t numcol, numrow;
	vector<vector<double>> fileGammaDaggervv;
	vector<double> fileFrequencyv;

	string file("/Users/drvdputt/GasModule/git/dat/t3_elec_reformat.ascii");
	ifstream input(file);
	if (!input) throw std::runtime_error("File " + file + " not found.");

	// The line number will not count any lines starting with #
	int lineNr = 0;
	string line;
	while(getline(input, line))
	{
		istringstream iss(line);
		if(iss.peek() != '#')
		{
			if (lineNr == 0)
			{
				iss >> numcol >> numrow;
				fileGammaDaggervv.resize(numrow, std::vector<double>(numcol, 0.));
			}
			if (lineNr == 1)
			{
				double log10T;
				while (iss >> log10T)
				{
					_logTemperaturev.push_back(log10T);
				}
				cout << endl;
			}
			if (lineNr > 1)
			{
				int flag;
				double energy;
				iss >> flag >> energy;

				double frequency = energy * Constant::RYDBERG / Constant::PLANCK;
				fileFrequencyv.push_back(frequency);
				if (flag) _thresholdv.push_back(frequency);

				auto it = fileGammaDaggervv[lineNr - 2].begin();
				double gammaDagger;
				while (iss >> gammaDagger)
				{
					*it = gammaDagger;
					it++;
				}
			}
			lineNr++;
		}
	}

	ofstream out;
	out.open("/Users/drvdputt/GasModule/bin/loadedContinuum.dat");
	for (auto& rowv : fileGammaDaggervv)
	{
		for(double d : rowv) out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// Now resample this data according to the wavelength grid
	_gammaDaggervv.resize(_wavelengthv.size(), vector<double>(numcol, 0));

	// First, from the wavelength grid, construct the desired frequency grid
	// (we may switch to frequencies for everything in the future, as these conversions are often needed)
	vector<double> frequencyv;
	frequencyv.reserve(_wavelengthv.size());
	for (auto rit = _wavelengthv.rbegin(); rit != _wavelengthv.rend(); rit++)
		frequencyv.push_back(Constant::LIGHT / *rit);

	cout << "frequency range from file: " << fileFrequencyv[0] << " to " << fileFrequencyv[fileFrequencyv.size()-1] << endl;
	cout << "frequency range: " << frequencyv[0] << " to " << frequencyv[frequencyv.size()-1] << endl;

	// Then, apply a linear interpolation on every column
	for (size_t col = 0; col < numcol; col++)
	{
		// Extract the column
		vector<double> column;
		column.reserve(numrow);
		for (size_t row = 0; row < numrow; row++) column.push_back(fileGammaDaggervv[row][col]);

		// Resample it
		const vector<double>& column_resampled = NumUtils::interpol<double>(column, fileFrequencyv, frequencyv, -1, -1);

		// And copy it over
		for (size_t row = 0; row < _gammaDaggervv.size(); row++)
			_gammaDaggervv[row][col] = column_resampled[row];
	}

	// DEBUG: print out the interpolated table
	out.open("/Users/drvdputt/GasModule/bin/interpolatedContinuum.dat");
	for (auto& rowv : _gammaDaggervv)
	{
		for(double d : rowv) out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// Now reverse the order to make the row indices correspond to wavelength indices
	reverse(_gammaDaggervv.begin(), _gammaDaggervv.end());
	// and convert the threshold frequencies to wavelengths
	for (double& t : _thresholdv) t = Constant::LIGHT / t;
	reverse(_thresholdv.begin(), _thresholdv.end());
}

void GasSpecies::solveBalance(double n, const vector<double>& isrf)
{
	_n = n;

	// Initial guess for the temperature
	_T = 500;

	calculateDensities(_T, isrf);
}

vector<double> GasSpecies::emissivity() const
{
	vector<double> result;
	const vector<double>& lineEmv = _levels.calculateEmission();
	const vector<double>& contEmv = continuumEmissivity(_T);
	for (size_t i = 0; i < _wavelengthv.size(); i++)
	{
		result.push_back(lineEmv[i] + contEmv[i]);
	}
	return result;
}

vector<double> GasSpecies::opacity() const
{
	vector<double> result;
	result.reserve(_wavelengthv.size());
	const vector<double>& lineOp = _levels.calculateOpacity();
	for (size_t i = 0; i < _wavelengthv.size(); i++)
	{
		double ionizOp_i = _n * (1 - _ionizedFraction) * Ionization::crossSection(_wavelengthv[i]);
		result.push_back(ionizOp_i + lineOp[i]);
	}
	return result;
}

void GasSpecies::calculateDensities(double T, const vector<double>& isrf)
{
	_ionizedFraction = Ionization::ionizedFraction(_n, T, _wavelengthv, isrf);

	cout << "Ionized fraction = " << _ionizedFraction << endl;

	double nAtomic = _n * (1 - _ionizedFraction);
	// For now we can sum over the collision partners, as TwoLevel threats all collision parnters in the same way
	// neutral + ion + electron densities
	double nTotal = (1 + _ionizedFraction) * _n;
	double np_ne_alpha = _n*_n*_ionizedFraction*_ionizedFraction * Ionization::recombinationRate(T);
	_levels.doLevels(nAtomic, nTotal, T, isrf, np_ne_alpha);
}

vector<double> GasSpecies::continuumEmissivity(double T) const
{
	double logT = log10(T);

	if (logT > _logTemperaturev.back() || logT < _logTemperaturev[0])
	{
		cout << "Warning: temperature " << T << " is outside of data range for free-bound continuum" << endl;
		return vector<double>(_wavelengthv.size(), 0.);
	}

	vector<double> result;
	result.reserve(_wavelengthv.size());

	// Find the grid point to the right of the requested log-temperature
	size_t iRight = NumUtils::index(logT, _logTemperaturev);
	// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
	double wRight = (T - _logTemperaturev[iRight - 1]) / (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

	// We will use equation (1) of Ercolano and Storey 2006 to remove the normalization of the data
	double Ttothe3_2 = pow(T, 1.5);
	double kT = Constant::BOLTZMAN * T;
	// "the nearest threshold of lower energy" i.e. then next threshold wavelength
	double tWav = 0;
	double tE;

	for (size_t iWav = 0; iWav < _wavelengthv.size(); iWav++)
	{
		double wav = _wavelengthv[iWav];

		// Interpolate gamma^dagger linearly in log T space
		double gammaDagger = (_gammaDaggervv[iWav][iRight - 1] * (1 - wRight)
				+ _gammaDaggervv[iWav][iRight] * wRight);

		// Skip over zero data, or when we are past the last threshold
		if (!gammaDagger || wav > _thresholdv.back())
		{
			result.push_back(0.);
		}
		else
		{
			// find a new tWav every time we pass a new threshold
			if (wav > tWav)
			{
				tWav = *upper_bound(_thresholdv.begin(), _thresholdv.end(), wav);
				tE = Constant::PLANCKLIGHT / tWav;
			}
			double E = Constant::PLANCKLIGHT / wav;
			double normalizationFactor = 1.e34 * Ttothe3_2 * exp((E - tE) / kT);
			result.push_back(gammaDagger / normalizationFactor);
		}
	}
	return result;
}
