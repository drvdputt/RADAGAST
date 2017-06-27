#include "FreeBound.h"
#include "Constants.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "SpecialFunctions.h"
#include "Table.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <cmath>
#include <fstream>

using namespace std;

FreeBound::FreeBound(const Array& frequencyv) : _frequencyv(frequencyv)
{
	/* Read the free-bound continuum data from Ercolano and Storey (2006). Adapted from NEBULAR
	   source code by M. Schirmer (2016). */
	vector<double> fileFrequencyv;
	vector<vector<double>> fileGammaDaggervv;
	vector<double> fileThresholdv;

	string file(repoRoot + "/dat/t3_elec_reformat.ascii");
	readData(file, fileFrequencyv, fileThresholdv, _logTemperaturev, fileGammaDaggervv);

	_thresholdv = Array(fileThresholdv.data(), fileThresholdv.size());

	size_t numcol = _logTemperaturev.size();
	size_t numrow = fileFrequencyv.size();

	// Now resample this data according to the frequency grid
	size_t numFreq = _frequencyv.size();
	_gammaDaggervv.resize(numFreq, numcol);

	DEBUG("frequency range from file: " << fileFrequencyv[0] << " to " << fileFrequencyv.back()
	                                    << endl);

	DEBUG("frequency range: " << _frequencyv[0] << " to " << _frequencyv[_frequencyv.size() - 1]
	                          << endl);

	/* Then, apply a linear interpolation across the frequencies (rows) for every temperature
	   (column) */
	for (size_t col = 0; col < numcol; col++)
	{
		// Extract the column
		vector<double> column;
		column.reserve(numrow);
		for (size_t row = 0; row < numrow; row++)
			column.push_back(fileGammaDaggervv[row][col]);

		// Resample it
		const vector<double>& column_resampled = NumUtils::interpol<double>(
		                column, fileFrequencyv,
		                vector<double>(begin(frequencyv), end(frequencyv)), -1, -1);

		// And copy it over
		for (size_t row = 0; row < numFreq; row++)
			_gammaDaggervv(row, col) = column_resampled[row];
	}

#ifdef DEBUG_CONTINUUM_DATA
	// DEBUG: print out the table as read from the file
	ofstream out = IOTools::ofstreamFile("freebound/loadedContinuum.dat");
	for (size_t iNu = 0; iNu < fileGammaDaggervv.size(); iNu++)
	{
		for (double d : fileGammaDaggervv[iNu])
			out << scientific << d << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: print out the interpolated table
	out = IOTools::ofstreamFile("freebound/interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (size_t iT = 0; iT < _gammaDaggervv.size(1); iT++)
			out << scientific << _gammaDaggervv(iNu, iT) << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: Test the temperature interpolation function (at least a copy pasta of a part)
	out = IOTools::ofstreamFile("freebound/bi-interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (double logT = 2; logT < 5; logT += 0.01)
		{
			// Find the grid point to the right of the requested log-temperature
			size_t iRight = NumUtils::index(logT, _logTemperaturev);
			/* The weight of the point to the right (= 1 if T is Tright, = 0 if T is
			   Tleft). */
			double wRight = (logT - _logTemperaturev[iRight - 1]) /
			                (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

			// Interpolate gamma^dagger linearly in log T space
			double gammaDagger = (_gammaDaggervv(iNu, iRight - 1) * (1 - wRight) +
			                      _gammaDaggervv(iNu, iRight) * wRight);
			out << scientific << gammaDagger << "\t";
		}
		out << endl;
	}
	out.close();
#endif /* DEBUG_CONTINUUM_DATA */

	DEBUG("Constructed FreeBound" << endl);
}

void FreeBound::readData(string file, vector<double>& fileFrequencyv,
                         vector<double>& fileThresholdv, vector<double>& fileTemperaturev,
                         vector<vector<double>>& fileGammaDaggervv) const
{
	ifstream input(IOTools::ifstreamFile(file));

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
	DEBUG("Successfully read freebound data in " << file << endl);
}

void FreeBound::addEmissionCoefficientv(double T, Array& gamma_nuv) const
{
	double logT = log10(T);

	if (logT > _logTemperaturev.back() || logT < _logTemperaturev[0])
	{
		DEBUG("Warning: temperature "
		      << T << "K is outside of data range for free-bound continuum" << endl);
		return;
	}

	// Find the grid point to the right of the requested log-temperature
	size_t iRight = NumUtils::index(logT, _logTemperaturev);
	size_t iLeft = iRight - 1;
	// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
	double wRight = (logT - _logTemperaturev[iLeft]) /
	                (_logTemperaturev[iRight] - _logTemperaturev[iLeft]);

	// We will use equation (1) of Ercolano and Storey 2006 to remove the normalization of the
	// data
	double Ttothe3_2 = pow(T, 3. / 2.);
	double kT = Constant::BOLTZMAN * T;
	// "the nearest threshold of lower energy"
	// i.e. the frequency in _thresholdv lower than the current one
	size_t iThreshold = 0;
	double tE = 0;
#ifdef DEBUG_CONTINUUM_DATA
	// TODO: Notice that this writes out the result every time this function is called. This
	// functionality needs to be moved someplace else.
	ofstream out = IOTools::ofstreamFile("freebound/gammanufb.dat");
#endif

	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double freq = _frequencyv[iFreq];

		// Interpolate gamma^dagger linearly in log T space
		double gammaDagger = _gammaDaggervv(iFreq, iLeft) +
		                     wRight * (_gammaDaggervv(iFreq, iRight) -
		                               _gammaDaggervv(iFreq, iLeft));

		double lastThresh = _thresholdv[_thresholdv.size() - 1];
		// Skip over zero data, or when we are below the first threshold
		if (!gammaDagger || freq < _thresholdv[0])
		{
		}
		else
		{
			// If the last threshold hasn't been passed yet, check if we have passed the
			// next (i + 1)
			if (freq < lastThresh)
			{
				// This block must only be executed when a new threshold is passed
				if (freq > _thresholdv[iThreshold + 1])
				{
					// find the next threshold of lower frequency
					// (don't just pick the next one, as this wouldn't work with
					// very coarse grids)
					iThreshold = TemplatedUtils::index(freq, _thresholdv) - 1;
					tE = Constant::PLANCK * _thresholdv[iThreshold];
				}
			}
			// When we have just moved past the last threshold, set the index one last
			// time
			else if (iThreshold < _thresholdv.size() - 1)
			{
				iThreshold = _thresholdv.size() - 1;
				tE = Constant::PLANCK * _thresholdv[iThreshold];
			}
			double E = Constant::PLANCK * freq;

			double normalizationFactor = 1.e34 * Ttothe3_2 * exp((E - tE) / kT);
			double gammaNu = gammaDagger / normalizationFactor;
#ifdef DEBUG_CONTINUUM_DATA
			out << freq << "\t" << gammaNu / 1.e-40 << endl;
#endif
			gamma_nuv[iFreq] += gammaNu;
		}
		// Also add the ionizing freebound continuum, which apparently is not included in
		// the data used here This is easily calculated using equation 4.21 from Osterbrock
		// and the Milne relation, since we have an expression for the photoionization
		// cross-section anyways
		if (freq > Ionization::THRESHOLD)
		{
			double h2 = Constant::PLANCK * Constant::PLANCK;
			double nu3 = freq * freq * freq;
			double m3 = Constant::ELECTRONMASS * Constant::ELECTRONMASS *
			            Constant::ELECTRONMASS;
			double c2 = Constant::LIGHT * Constant::LIGHT;
			double u_nu = sqrt(2 / Constant::ELECTRONMASS * Constant::PLANCK *
			                   (freq - Ionization::THRESHOLD));
			double u2 = u_nu * u_nu;
			double a_nu = Ionization::crossSection(freq);
			double f_u_nu = SpecialFunctions::maxwellBoltzman(u_nu, T,
			                                                  Constant::ELECTRONMASS);
			double ionizingContinuum = h2 * h2 * nu3 / u2 / m3 / c2 * a_nu * f_u_nu;
			gamma_nuv[iFreq] += ionizingContinuum;
		}
	}
#ifdef DEBUG_CONTINUUM_DATA
	out.close();
#endif
}
