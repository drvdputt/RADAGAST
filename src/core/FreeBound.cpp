#include "FreeBound.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "Testing.h"

#include <cmath>
#include <fstream>

using namespace std;

FreeBound::FreeBound()
{
	/* Read the free-bound continuum data from Ercolano and Storey (2006). Adapted from
	   NEBULAR source code by M. Schirmer (2016). */

	// Use vectors first, for flexibility
	vector<double> fileFrequencyv;
	vector<vector<double>> fileGammaDaggervv;
	vector<double> fileThresholdv;
	vector<double> fileTemperaturev;

	readData("dat/t3_elec_reformat.ascii", fileFrequencyv, fileThresholdv, fileTemperaturev,
	         fileGammaDaggervv);

	// Switch to Array here, because that's what we use everywhere in the code
	_frequencyv = Array(fileFrequencyv.data(), fileFrequencyv.size());
	_thresholdv = Array(fileThresholdv.data(), fileThresholdv.size());
	_logTemperaturev = Array(fileTemperaturev.data(), fileTemperaturev.size());

	// The number of temperatures
	size_t numcol = _logTemperaturev.size();
	// The number of frequencies
	size_t numrow = _frequencyv.size();

	// Copy everything into this Table for speed
	_gammaDaggervv.resize(numrow, numcol);
	for (size_t row = 0; row < numrow; row++)
	{
		double* table_row_begin = &_gammaDaggervv(row, 0);
		copy(fileGammaDaggervv[row].begin(), fileGammaDaggervv[row].end(),
		     table_row_begin);
	}

	DEBUG("Constructed FreeBound" << endl);
// #define DEBUG_CONTINUUM_DATA
#ifdef DEBUG_CONTINUUM_DATA
	// DEBUG: print out the table
	ofstream out = IOTools::ofstreamFile("freebound/loadedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (size_t iT = 0; iT < _gammaDaggervv.size(1); iT++)
			out << scientific << _gammaDaggervv(iNu, iT) << '\t';
		out << endl;
	}
	out.close();

	// DEBUG: Test the temperature interpolation function (at least a copy pasta of a part)
	out = IOTools::ofstreamFile("freebound/T-interpolatedContinuum.dat");
	for (size_t iNu = 0; iNu < _gammaDaggervv.size(0); iNu++)
	{
		for (double logT = 2; logT < 5; logT += 0.01)
		{
			// Find the grid point to the right of the requested log-temperature
			size_t iRight = TemplatedUtils::index(logT, _logTemperaturev);
			/* The weight of the point to the right (= 1 if T is Tright, = 0 if T is
			   Tleft). */
			double wRight = (logT - _logTemperaturev[iRight - 1]) /
			                (_logTemperaturev[iRight] -
			                 _logTemperaturev[iRight - 1]);

			// Interpolate gamma^dagger linearly in log T space
			double gammaDagger = (_gammaDaggervv(iNu, iRight - 1) * (1 - wRight) +
			                      _gammaDaggervv(iNu, iRight) * wRight);
			out << scientific << gammaDagger << "\t";
		}
		out << endl;
	}
	out.close();

	// DEBUG: test the gammanu function to obtain a figure like in the nebular paper
	Array testFrequencyv = Testing::generateGeometricGridv(1000, _frequencyv[0], 1e16);

	Array gammaNuv(testFrequencyv.size());
	addEmissionCoefficientv(10000., testFrequencyv, gammaNuv);

	out = IOTools::ofstreamFile("freebound/gammanufb.dat");
	for (size_t iNu = 0; iNu < testFrequencyv.size(); iNu++)
		out << testFrequencyv[iNu] << "\t" << gammaNuv[iNu] / 1.e-40 << endl;
	out.close();
#endif /* DEBUG_CONTINUUM_DATA */
}

void FreeBound::readData(string file, vector<double>& fileFrequencyv,
                         vector<double>& fileThresholdv, vector<double>& fileTemperaturev,
                         vector<vector<double>>& fileGammaDaggervv) const
{
	ifstream input(IOTools::ifstreamRepoFile(file));

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
				fileGammaDaggervv.resize(numrow,
				                         std::vector<double>(numcol, 0.));
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

				double frequency =
				                energy * Constant::RYDBERG / Constant::PLANCK;
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

void FreeBound::addEmissionCoefficientv(double T, const Array& eFrequencyv,
                                        Array& gamma_nuv) const
{
	double logT = log10(T);

	if (logT > _logTemperaturev[_logTemperaturev.size() - 1] || logT < _logTemperaturev[0])
	{
		DEBUG("Warning: temperature "
		      << T << "K is outside of data range for free-bound continuum" << endl);
		return;
	}

	// Find the grid point to the right of the requested log-temperature
	size_t iRight = TemplatedUtils::index(logT, _logTemperaturev);
	size_t iLeft = iRight - 1;
	// The weight of the point to the right (= 1 if T is Tright, = 0 if T is Tleft)
	double wRight = (logT - _logTemperaturev[iLeft]) /
	                (_logTemperaturev[iRight] - _logTemperaturev[iLeft]);

	// Interpolate the data for this temperature
	Array t_interpolated_gammaDaggerv(_frequencyv.size());
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		// If this loop is (to a significant degree) slow, consider transposing the
		// gammaDaggervv table, to do this loop contiguously

		// Interpolate gamma^dagger linearly in log T space
		t_interpolated_gammaDaggerv[iFreq] = _gammaDaggervv(iFreq, iLeft) +
		                                     wRight * (_gammaDaggervv(iFreq, iRight) -
		                                               _gammaDaggervv(iFreq, iLeft));
	}

	// Put onto the right frequency grid
	Spectrum gammaDaggerSp(_frequencyv, t_interpolated_gammaDaggerv);
	// Think about whether we should pick 'evaluate' or 'binned' here. (if a bin spans over
	// a of the normalization threshold, we run into trouble...).
	Array nu_interpolated_gammaDaggerv = gammaDaggerSp.binned(eFrequencyv);

	/* We will use equation (1) of Ercolano and Storey 2006 to remove the normalization of
	   the data. */
	double Ttothe3_2 = pow(T, 3. / 2.);
	double kT = Constant::BOLTZMAN * T;
	/* "the nearest threshold of lower energy" i.e. the frequency in _thresholdv lower than
	   the current one. */
	size_t iThreshold = 0;
	double lastThresh = _thresholdv[_thresholdv.size() - 1];
	double tE = 0;

	// Only now apply the exponential factors (important for interpolation quality)
	for (size_t i = 0; i < eFrequencyv.size(); i++)
	{
		double freq = eFrequencyv[i];

		// Skip over zero data, or when we are below the first threshold
		if (nu_interpolated_gammaDaggerv[i] && freq >= _thresholdv[0])
		{
			// If the last threshold hasn't been passed yet, check if we have passed
			// the next (i + 1)
			if (freq < lastThresh)
			{
				// This block must only be executed when a new threshold is
				// passed
				if (freq > _thresholdv[iThreshold + 1])
				{
					// find the next threshold of lower frequency (don't
					// just pick the next one, as this wouldn't work with
					// very coarse grids)
					iThreshold = TemplatedUtils::index(freq, _thresholdv) -
					             1;
					tE = Constant::PLANCK * _thresholdv[iThreshold];
				}
			}
			// When we have just moved past the last threshold, set the index one
			// last time
			else if (iThreshold < _thresholdv.size() - 1)
			{
				iThreshold = _thresholdv.size() - 1;
				tE = Constant::PLANCK * _thresholdv[iThreshold];
			}
			double E = Constant::PLANCK * freq;

			double normalizationFactor = 1.e34 * Ttothe3_2 * exp((E - tE) / kT);
			nu_interpolated_gammaDaggerv[i] /= normalizationFactor;
		}
	}

	gamma_nuv += nu_interpolated_gammaDaggerv;

	// Also add the ionizing freebound continuum, which apparently is not included in
	// the data used here This is easily calculated using equation 4.21 from Osterbrock
	// and the Milne relation, since we have an expression for the photoionization
	// cross-section anyways
	for (size_t iFreq = 0; iFreq < eFrequencyv.size(); iFreq++)
	{
		double freq = eFrequencyv[iFreq];
		if (freq > Ionization::THRESHOLD)
			gamma_nuv[iFreq] += ionizingContinuum(T, freq);
	}
}

double FreeBound::ionizingContinuum(double T, double frequency)
{
	double h2 = Constant::PLANCK * Constant::PLANCK;
	double nu3 = frequency * frequency * frequency;
	double m3 = Constant::ELECTRONMASS * Constant::ELECTRONMASS * Constant::ELECTRONMASS;
	double c2 = Constant::LIGHT * Constant::LIGHT;
	double u_nu = sqrt(2 / Constant::ELECTRONMASS * Constant::PLANCK *
	                   (frequency - Ionization::THRESHOLD));
	double u2 = u_nu * u_nu;
	double a_nu = Ionization::crossSection(frequency);
	double f_u_nu = SpecialFunctions::maxwellBoltzman(u_nu, T, Constant::ELECTRONMASS);
	return h2 * h2 * nu3 / u2 / m3 / c2 * a_nu * f_u_nu;
}
