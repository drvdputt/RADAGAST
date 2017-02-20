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
	: _wavelengthv(wavelengthv), _n(0), _ionizedFraction(0), _levels(wavelengthv)
{
	// Read the free-bound continuum data from Ercolano and Storey (2006)
	// Adapted from NEBULAR source code by M. Schirmer (2016)
	int numcol, numrow;
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
				fileFrequencyv.push_back(energy * Constant::RYDBERG / Constant::PLANCK);

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

	// First, construct the corresponding frequency grid
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
}

void GasSpecies::solveBalance(double n, const vector<double>& isrf)
{
	_n = n;

	// Initial guess for the temperature
	double T = 1000;

	calculateDensities(T, isrf);
}

vector<double> GasSpecies::emissivity() const
{
	return _levels.calculateEmission();
}

vector<double> GasSpecies::opacity() const
{
	return _levels.calculateOpacity();
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

vector<double> GasSpecies::continuumEmissivity(double T)
{
	// Interpolate linearly between the log-temperature grid points
	double logT = log10(T);
	size_t iRight = NumUtils::index(logT, _logTemperaturev);
	bool edge = false;
	if (iTleft == 0)
	{
		edge = true;
	}
	if (iTleft == _logTemperaturev.size() - 1)
	{
		edge = true;
	}
	double Tleft = _logTemperaturev[iTleft];
	double Tright =

	vector<double> result(_wavelengthv.size());
	for (size_t iWav = 0; iWav < _wavelengthv.size(); iWav++)
	{
		const vector<double>& gammaDaggerv = _gammaDaggervv[iWav];


	}
	return result
}
