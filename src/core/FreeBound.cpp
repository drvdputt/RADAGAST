#include "FreeBound.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include <gsl/gsl_pow_int.h>
#include <cmath>
#include <fstream>

using namespace std;

namespace GasModule
{
    FreeBound::FreeBound()
    {
        /* Read the free-bound continuum data from Ercolano and Storey (2006). Adapted from
	   NEBULAR source code by M. Schirmer (2016). */

        // Use vectors first, for flexibility
        vector<double> fileFrequencyv;
        vector<vector<double>> fileGammaDaggervv;
        vector<double> fileThresholdv;
        vector<double> fileTemperaturev;

        readData("dat/t3_elec_reformat.ascii", fileFrequencyv, fileThresholdv, fileTemperaturev, fileGammaDaggervv);

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
            copy(fileGammaDaggervv[row].begin(), fileGammaDaggervv[row].end(), table_row_begin);
        }

        DEBUG("Constructed FreeBound" << endl);
        if (Options::freebound_debugData)
        {
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
                    double wRight = (logT - _logTemperaturev[iRight - 1])
                                    / (_logTemperaturev[iRight] - _logTemperaturev[iRight - 1]);

                    // Interpolate gamma^dagger linearly in log T space
                    double gammaDagger =
                        (_gammaDaggervv(iNu, iRight - 1) * (1 - wRight) + _gammaDaggervv(iNu, iRight) * wRight);
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
        }
    }

    void FreeBound::readData(string file, vector<double>& fileFrequencyv, vector<double>& fileThresholdv,
                             vector<double>& fileTemperaturev, vector<vector<double>>& fileGammaDaggervv) const
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
                    fileGammaDaggervv.resize(numrow, std::vector<double>(numcol, 0.));
                }
                if (lineNr == 1)
                {
                    double log10T;
                    while (iss >> log10T) fileTemperaturev.push_back(log10T);
                }
                if (lineNr > 1)
                {
                    int flag;
                    double energy;
                    iss >> flag >> energy;

                    double frequency = energy * Constant::RYDBERG / Constant::PLANCK;
                    fileFrequencyv.push_back(frequency);
                    if (flag) fileThresholdv.push_back(frequency);

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

    void FreeBound::addEmissionCoefficientv(double T, const Array& eFrequencyv, Array& gamma_nuv) const
    {
        double logT = log10(T);
        if (logT > _logTemperaturev[_logTemperaturev.size() - 1] || logT < _logTemperaturev[0])
        {
            DEBUG("Warning: temperature " << T << "K is outside of data range for free-bound continuum" << endl);
            return;
        }

        // The data are rather coarse (only ~ 50 frequency points). According to the paper, we
        // should interpolate linearly it in both temperature and frequency space before
        // applying the normalization (= exponentiating).
        size_t tLeft, tRight;
        std::tie(tLeft, tRight) = TemplatedUtils::saneIndexPair(logT, _logTemperaturev);

        // Interpolate the data between these two temperature points, for all frequencies
        int numOutputFreqs = eFrequencyv.size();
        if (numOutputFreqs < 3) Error::runtime("eFrequencyv should have at least 3 points");
        Array interpolated_gammaDaggerv(numOutputFreqs);

        // Attention: do not integrate over bins here. It causes problems when a bin overlaps
        // one of the normalization thresholds. Instead, we will just evaluate the spectrum at
        // the frequency points using linear interpolation.

        // Start with data points [0, 1] and update when freq goes above current upper point
        size_t nuLo{0}, nuHi{1};
        double minDataFreq = _frequencyv[0];
        double maxDataFreq = _frequencyv[_frequencyv.size() - 1];
        for (size_t i = 0; i < numOutputFreqs; i++)
        {
            double nu = eFrequencyv[i];
            if (nu < minDataFreq || maxDataFreq < nu)
            {
                interpolated_gammaDaggerv[i] = 0;
                continue;
            }

            // Find new nuLo and nuHi pair when necessary
            if (nu > _frequencyv[nuHi]) std::tie(nuLo, nuHi) = TemplatedUtils::saneIndexPair(nu, _frequencyv);

            // Interpolate gamma^dagger bilinearly in logT-nu space.
            interpolated_gammaDaggerv[i] = TemplatedUtils::interpolateBilinear(
                logT, nu, _logTemperaturev[tLeft], _logTemperaturev[tRight], _frequencyv[nuLo], _frequencyv[nuHi],
                _gammaDaggervv(nuLo, tLeft), _gammaDaggervv(nuLo, tRight), _gammaDaggervv(nuHi, tLeft),
                _gammaDaggervv(nuHi, tRight));
        }

        // Use equation (1) of Ercolano and Storey 2006 to remove the normalization of the data.
        double Ttothe3_2 = pow(T, 3. / 2.);
        double kT = Constant::BOLTZMAN * T;
        // Keep track of "the nearest threshold of lower energy" i.e. the frequency in
        // _thresholdv lower than the current one.
        int iThreshold = 0;
        int iLastThresh = _thresholdv.size() - 1;
        double lastThresh = _thresholdv[iLastThresh];
        double tE = 0;

        // Apply the exponential factors
        for (size_t i = 0; i < numOutputFreqs; i++)
        {
            double nu = eFrequencyv[i];

            // Skip over zero data, or when we are below the first threshold
            if (interpolated_gammaDaggerv[i] <= 0 || nu < _thresholdv[0]) continue;

            // If the last threshold hasn't been passed yet, we need to check if we have passed
            // the next one and adjust the normalization parameters
            if (nu < lastThresh)
            {
                if (nu > _thresholdv[iThreshold + 1])
                {
                    // find the next threshold of lower frequency (don't just pick the next one, as
                    // this wouldn't work with very coarse grids)
                    iThreshold = TemplatedUtils::index(nu, _thresholdv) - 1;
                }
                // else do nothing
            }
            // When we are past the last threshold, use the energy of the last one
            else
                iThreshold = iLastThresh;

            tE = Constant::PLANCK * _thresholdv[iThreshold];
            double deltaE = Constant::PLANCK * nu - tE;
            double normalizationFactor = 1.e34 * Ttothe3_2 * exp(deltaE / kT);
            interpolated_gammaDaggerv[i] /= normalizationFactor;
        }

        // Add the interpolated data to the total SED
        gamma_nuv += interpolated_gammaDaggerv;

        // Also add the ionizing freebound continuum, which apparently is not included in the
        // data used here. This is easily calculated using equation 4.21 from Osterbrock and the
        // Milne relation, since we have an expression for the photoionization cross-section
        // anyways
        for (size_t i = 0; i < numOutputFreqs; i++) gamma_nuv[i] += ionizingContinuumCoefficient(T, eFrequencyv[i]);
    }

    double FreeBound::ionizingContinuumCoefficient(double T, double frequency)
    {
        if (frequency < Ionization::THRESHOLD) return 0.;

        // I found it non-trivial to find this whole formula written out somewhere. This
        // conference paper has it: R.R. Sheeba et al 2019 J. Phys.: Conf. Ser.1289 012041. To
        // derive it yourself: start from power_density_nu dnu = hnu * np * ne * ve * sigma(ve)
        // * f(ve) dve. Every electron speed ve corresponds just 1 hnu = me ve^2 / 2 + I
        // (kinetic + ionization energy). Apply the milne relation to replace sigma(ve) (the
        // capture cross section), and fill in the maxwell distribution f(ve). Then sustitute ve
        // = sqrt(2 (hnu - I) / me). Note that the weight ratio in the milne relation is 2; H1s
        // has double the multiplicity of a proton, because of the electron spin.
        const double factors = 4 * sqrt(2. / Constant::PI) * gsl_pow_4(Constant::PLANCK) / Constant::LIGHT2;
        double kT = Constant::BOLTZMAN * T;
        double hnuDelta = Constant::PLANCK * (frequency - Ionization::THRESHOLD);
        double emCoeff_nu = factors * pow(Constant::ELECTRONMASS * kT, -1.5) * exp(-hnuDelta / kT)
                            * Ionization::crossSection(frequency) * gsl_pow_3(frequency);
        return emCoeff_nu;  // is multiplied with ne * ne / 4pi elsewhere in the code
    }
}
