#include "LevelCoefficients.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Functions.hpp"
#include "Options.hpp"
#include "TemplatedUtils.hpp"

using namespace std;

namespace GasModule
{
    LevelCoefficients::LevelCoefficients(double mass) : _mass(mass)
    {
        // // Do a sanity check: All active transitions must be downward ones in energy
        // forActiveLinesDo([&](size_t upper, size_t lower) {
        // 	if (_ev(upper) - _ev(lower) < 0)
        // 	{
        // 		cout << "Eu - El = " << _ev(upper) - _ev(lower) << endl;
        // 		cout << "u l = " << upper << " " << lower
        // 		     << " Aul = " << _avv(upper, lower) << endl;
        // 		Error::runtime("There is an upward A-coefficient. This can't be "
        // 		               "correct");
        // 	}
        // });
    }

    void LevelCoefficients::setConstants(const EVector& ev, const EVector& gv, const EMatrix& avv,
                                         const EMatrix& extraAvv)
    {
        _numLv = ev.size();
        _ev = ev;
        _gv = gv;
        _avv = avv;
        _extraAvv = extraAvv;
    }

    void LevelCoefficients::lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const
    {
        // Count the number of active lines
        numLines = 0;
        forActiveLinesDo([&](size_t, size_t) { numLines++; });
        lineFreqv.resize(numLines);
        naturalLineWidthv.resize(numLines);

        // Calculate the natural line width for these transitions
        int index = 0;
        forActiveLinesDo([&](size_t upper, size_t lower) {
            lineFreqv[index] = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
            naturalLineWidthv[index] = (_avv(upper, lower) + _extraAvv(upper, lower)) / Constant::FPI;
            index++;
        });
    }

    EMatrix LevelCoefficients::prepareAbsorptionMatrix(const Spectrum& meanIntensity, double T,
                                                       const EMatrix& Cvv) const
    {
        EMatrix BPvv = EMatrix::Zero(_numLv, _numLv);
        double spectrumMax = meanIntensity.valMax();
        forActiveLinesDo([&](size_t upper, size_t lower) {
            // Calculate Pij for the lower triangle (= stimulated emission)
            LineProfile lp = lineProfile(upper, lower, T, Cvv);
            double linePower = lp.integrateSpectrum(meanIntensity, spectrumMax);

            if (Options::levelcoefficients_reportSpecIntegral)
            {
                // Compare above integral to explicit calculation using the whole wavelength
                // range with many points. Useful check in case of steep slopes in the spectrum,
                // or very wide wings of the line.
                auto f = [&](double x) -> double { return lp(x) * meanIntensity.evaluate(x); };

                size_t many_points = 1e6;
                double manualIntegral = TemplatedUtils::integrateFunction<double>(f, meanIntensity.freqMin(),
                                                                                  meanIntensity.freqMax(), many_points);

                double ratio = manualIntegral / linePower;

                if (abs(ratio - 1.) > 1.e-6) cout << linePower << "\t MR:" << ratio << endl;
            }

            // Multiply by Bij in terms of Aij, valid for i > j
            double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
            BPvv(upper, lower) = Constant::CSQUARE_TWOPLANCK / (nu_ij * nu_ij * nu_ij) * _avv(upper, lower) * linePower;

            /* Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij = Pji. */
            BPvv(lower, upper) = _gv(upper) / _gv(lower) * BPvv(upper, lower);
        });
        return BPvv;
    }

    EVector LevelCoefficients::solveBoltzmanEquations(double T) const
    {
        double eMin = _ev.minCoeff();
        // Degeneracy factor
        EVector pv{_gv};
        // Partition function
        double pSum{0};
        double kT = Constant::BOLTZMAN * T;
        for (int i = 0; i < _ev.size(); i++)
        {
            // Exponential factor
            pv(i) *= exp((eMin - _ev(i)) / kT);
            // Normalization
            pSum += pv(i);
        }
        return pv / pSum;
    }

    void LevelCoefficients::forActiveLinesDo(function<void(size_t ini, size_t fin)> thing) const
    {
        // Execute the same function for all transitions that are optically active.
        for (size_t fin = 0; fin < _numLv; fin++)
            for (size_t ini = 0; ini < _numLv; ini++)
                if (_avv(ini, fin)) thing(ini, fin);
    }

    double LevelCoefficients::lineIntensityFactor(size_t upper, size_t lower, double nu) const
    {
        return (_ev(upper) - _ev(lower)) / Constant::FPI * nu * _avv(upper, lower);
    }

    double LevelCoefficients::lineOpacityFactor(size_t upper, size_t lower, double nu, double nl) const
    {
        double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
        double constantFactor =
            Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij / nu_ij * _avv(upper, lower);
        double densityFactor = nl * _gv(upper) / _gv(lower) - nu;
        double result = constantFactor * densityFactor;
        return result;
    }

    LineProfile LevelCoefficients::lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const
    {
        double nu0 = (_ev(upper) - _ev(lower)) / Constant::PLANCK;

        double decayRate = _avv(upper, lower) + _extraAvv(upper, lower) + Cvv(upper, lower)  // decay rate of top level
                           + Cvv(lower, upper);  // decay rate of bottom level
        // (stimulated emission doesn't count, as it causes no broadening)

        // Half the FWHM of the Lorentz
        double halfWidth = decayRate / Constant::FPI;

        // The standard deviation in frequency units. It is about half of the FWHM for a
        // Gaussian
        double sigma_nu = nu0 * Functions::thermalVelocityWidth(T, _mass) / Constant::LIGHT;

        return LineProfile(nu0, sigma_nu, halfWidth);
    }
}
