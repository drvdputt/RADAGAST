#include "FreeBound.h"
#include "flags.h"
#include "Constants.h"
#include "NumUtils.h"
#include "Table.h"

#include <cmath>
#include <fstream>

using namespace std;

FreeBound::FreeBound(const vector<double>& frequencyv) :
		_frequencyv(frequencyv)
{
// Read the free-bound continuum data from Ercolano and Storey (2006)
// Adapted from NEBULAR source code by M. Schirmer (2016)
	vector<double> fileFrequencyv;
	vector<vector<double>> fileGammaDaggervv;

	string file("/Users/drvdputt/GasModule/git/dat/t3_elec_reformat.ascii");
	readData(file, fileFrequencyv, _thresholdv, _logTemperaturev, fileGammaDaggervv);

	size_t numcol = _logTemperaturev.size();
	size_t numrow = fileFrequencyv.size();

	// Now resample this data according to the frequency grid
	size_t numFreq = _frequencyv.size();
	_gammaDaggervv.resize(numFreq, numcol);

	DEBUG("frequency range from file: " << fileFrequencyv[0] << " to " << fileFrequencyv.back() << endl);

	DEBUG("frequency range: " << _frequencyv[0] << " to " << _frequencyv.back() << endl);

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
		for (size_t row = 0; row < numFreq; row++)
			_gammaDaggervv(row, col) = column_resampled[row];
	}

#ifdef PRINT_CONTINUUM_DATA
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
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (size_t iT = 0; iT < _gammaDaggervv.size(1); iT++)
		out << scientific << _gammaDaggervv(iNu, iT) << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: Test the temperature interpolation function (at least a copy pasta of a part)
	out.open("/Users/drvdputt/GasModule/run/bi-interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (double logT = 2; logT < 5; logT += 0.01)
		{
			// Find the grid point to the right of the requested log-temperature
			size_t iRight = NumUtils::index(logT, _logTemperaturev);
			// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
			double wRight = (logT - _logTemperaturev[iRight - 1])
			/ (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

			// Interpolate gamma^dagger linearly in log T space
			double gammaDagger = (_gammaDaggervv(iNu, iRight - 1) * (1 - wRight)
					+ _gammaDaggervv(iNu, iRight) * wRight);
			out << scientific << gammaDagger << "\t";
		}
		out << endl;
	}
	out.close();
#endif /*PRINT_CONTINUUM_DATA*/
}

void FreeBound::readData(string file, vector<double>& fileFrequencyv, vector<double>& fileThresholdv,
		vector<double>& fileTemperaturev, vector<vector<double>>& fileGammaDaggervv) const
{
	ifstream input(file);
	if (!input)
		throw std::runtime_error("File " + file + " not found.");

	size_t numcol, numrow;

	// The line number will not count any lines starting with #
	int lineNr = 0;
	string line;
	while (getline(input, line))
	{
		istringstream iss(line);
		if (iss.peek() != '#')
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
					fileTemperaturev.push_back(log10T);
			}
			if (lineNr > 1)
			{
				int flag;
				double energy;
				iss >> flag >> energy;

				double frequency = energy * Constant::RYDBERG / Constant::PLANCK;
				fileFrequencyv.push_back(frequency);
				if (flag)
					fileThresholdv.push_back(frequency);

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
}

vector<double> FreeBound::emissionCoefficientv(double T) const
{
	double logT = log10(T);

	if (logT > _logTemperaturev.back() || logT < _logTemperaturev[0])
	{
#ifdef VERBOSE
		DEBUG(
				"Warning: temperature " << T << "K is outside of data range for free-bound continuum" << endl);
#endif
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
#ifdef PRINT_CONTINUUM_DATA
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/gammanu.dat");
#endif

	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];

		// Interpolate gamma^dagger linearly in log T space
		double gammaDagger = _gammaDaggervv(iFreq, iRight - 1) * (1 - wRight)
				+ _gammaDaggervv(iFreq, iRight) * wRight;

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
#ifdef PRINT_CONTINUUM_DATA
			out << freq << "\t" << gammaNu / 1.e-40 << endl;
#endif
			result.push_back(gammaNu);
		}
	}
#ifdef PRINT_CONTINUUM_DATA
	out.close();
#endif
	return result;
}

