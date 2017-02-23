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
	for (size_t iNu = 0; iNu < fileGammaDaggervv.size(); iNu++)
	{
		for(double d : fileGammaDaggervv[iNu]) out << scientific << d << '\t';
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
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(); iNu++)
	{
		for(double d : _gammaDaggervv[iNu]) out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// Now reverse the order to make the row indices correspond to wavelength indices
	reverse(_gammaDaggervv.begin(), _gammaDaggervv.end());
	// and convert the threshold frequencies to wavelengths
	for (double& t : _thresholdv) t = Constant::LIGHT / t;
	reverse(_thresholdv.begin(), _thresholdv.end());

	// Test the temperature interpolation function (at least a copy pasta of a part)
	out.open("/Users/drvdputt/GasModule/bin/bi-interpolatedContinuum.dat");
	for (size_t iWav = 0; iWav < _gammaDaggervv.size(); iWav++)
	{
		for (double logT = 2; logT < 5; logT += 0.01)
		{
			// Find the grid point to the right of the requested log-temperature
			size_t iRight = NumUtils::index(logT, _logTemperaturev);
			// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
			double wRight = (logT - _logTemperaturev[iRight - 1]) / (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

			// Interpolate gamma^dagger linearly in log T space
			double gammaDagger = (_gammaDaggervv[iWav][iRight - 1] * (1 - wRight)
					+ _gammaDaggervv[iWav][iRight] * wRight);
			out << scientific << gammaDagger << "\t";
		}
		out << endl;
	}
	out.close();
}

void GasSpecies::solveBalance(double n, double Tinit, const vector<double>& isrf)
{
	_n = n;

	// Initial guess for the temperature
	_T = Tinit;

	calculateDensities(_T, isrf);
}

vector<double> GasSpecies::emissivity() const
{
	vector<double> result;
	const vector<double>& lineEmv = _levels.calculateEmission();
	const vector<double>& contEmv = continuumEmissionCoeff(_T);
	for (size_t i = 0; i < _wavelengthv.size(); i++)
	{
		result.push_back(lineEmv[i] + _n*_n*_ionizedFraction*_ionizedFraction / Constant::FPI * contEmv[i]);
	}
	cout << "line / cont = "
			<< NumUtils::integrate<double>(_wavelengthv, lineEmv) << " / "
			<< _n*_n*_ionizedFraction*_ionizedFraction / Constant::FPI * NumUtils::integrate<double>(_wavelengthv, contEmv) << endl;
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
	cout << "Calculating state for T = " << T << "K" << endl;

	_ionizedFraction = Ionization::ionizedFraction(_n, T, _wavelengthv, isrf);

	cout << "Ionized fraction = " << _ionizedFraction << endl;

	double nAtomic = _n * (1 - _ionizedFraction);
	// For now we can sum over the collision partners, as TwoLevel threats all collision parnters in the same way
	// neutral + ion + electron densities
	double nTotal = (1 + _ionizedFraction) * _n;

	double ne_np = _n*_n*_ionizedFraction*_ionizedFraction;
	double alphaTotal = Ionization::recombinationRate(T);

	// approximations from Draine's book, p 138, valid for 3000 to 30000 K
	double T4 = T / 1.e4;
	double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4)); // yes, this is natural log
	double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
	cout << "alphaGround " << alphaGround << " alpha2p " << alpha2p << endl;

	vector<double> sourcev = {ne_np * alphaGround, ne_np * alpha2p};

	// The ionization rate calculation makes no distinction between the levels.
	// Therefore, the sink term is the same for each level. Moreover, total source = total sink
	// so we want sink*n0 + sink*n1 = source => sink = source / n because n0/n + n1/n = 1
	double sink = ne_np * alphaTotal / nAtomic;
	vector<double> sinkv = {sink, sink};

	_levels.doLevels(nAtomic, nTotal, T, isrf, sourcev, sinkv);
}

vector<double> GasSpecies::continuumEmissionCoeff(double T) const
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
	double wRight = (logT - _logTemperaturev[iRight - 1]) / (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

	// We will use equation (1) of Ercolano and Storey 2006 to remove the normalization of the data
	double Ttothe3_2 = pow(T, 3./2.);
	double kT = Constant::BOLTZMAN * T;
	// "the nearest threshold of lower energy" i.e. then next threshold wavelength
	double tWav = 0;
	double tE = _thresholdv[0];

	ofstream out;
	out.open("/Users/drvdputt/GasModule/bin/gammanu.dat");

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
				// (don't just pick the next one, as this wouldn't work with very coarse grids)
				tWav = *upper_bound(_thresholdv.begin(), _thresholdv.end(), wav);
				tE = Constant::PLANCKLIGHT / tWav;
			}
			double E = Constant::PLANCKLIGHT / wav;
			double normalizationFactor = 1.e34 * Ttothe3_2 * exp((E - tE) / kT);

			double gammaNu = gammaDagger / normalizationFactor;
			out << Constant::LIGHT/wav << "\t" << gammaNu / 1.e-40 << endl;

			// The order is already correct, but still need to convert to emissivity per wavelength interval
			result.push_back(gammaNu * Constant::LIGHT / wav/wav);
		}
	}
	out.close();
	return result;
}
