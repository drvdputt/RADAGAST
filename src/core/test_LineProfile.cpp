#include "Array.hpp"
#include "DoctestUtils.hpp"
#include "Functions.hpp"
#include "HFromFiles.hpp"
#include "LevelCoefficients.hpp"
#include "LineProfile.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "doctest.h"

using namespace RADAGAST;

namespace
{
    double _c = 4.;

    LineProfile lpFactory(const std::string& testCase)
    {
        if (testCase == "mixed")
            return LineProfile(_c, 0.1, 0.1);
        else if (testCase == "zero")
            return LineProfile(_c, 0, .1);
        else  // (testCase == "thermal")
            return LineProfile(_c, .2, .01);
    }

    Spectrum testSpectrum(double base)
    {
        size_t numFreq = 100;
        Array frequencyv(numFreq);
        double fmin = .01;
        double fmax = 10.;
        double step = (fmax - fmin) / numFreq;
        for (size_t i = 0; i < numFreq; i++) frequencyv[i] = fmin + i * step;
        Array spectrumv(base, numFreq);
        return Spectrum(frequencyv, spectrumv);
    }
}  // namespace

TEST_CASE("Test line profile integration over spectrum")
{
    // Generate a flat spectrum of value base
    double base = 5.;
    auto s = testSpectrum(base);

    std::string testCase;
    double e = 1.e-2;
    bool warn = false;

    SUBCASE("mixed guassian/lorentzian line")
    {
        testCase = "mixed";
        // not expected to work perfectly at the moment
        warn = true;
    }

    SUBCASE("mostly gaussian line") { testCase = "thermal"; }

    SUBCASE("zero temp") { testCase = "zero"; }

    auto lp = lpFactory(testCase);
    double integral = lp.integrateSpectrum(s, 0);
    // This integral should be about equal to base, if gridpoints are chosen well
    DoctestUtils::checkTolerance("test value", integral, base, e, warn);

    // we can test the binned addition here too
    Array freqv = s.frequencyv();
    Array flux = s.valuev();
    // last fector in addToSpectrum is integrated line flux
    double lineStrength = 2;
    double before = TemplatedUtils::integrate<double, Array, Array>(freqv, flux);
    lp.addToBinned(freqv, flux, lineStrength);
    double after = TemplatedUtils::integrate<double, Array, Array>(freqv, flux);

    // check if int(new flux) dnu = int(old flux) dnu + lineStrength
    DoctestUtils::checkTolerance("old flux + line", after, before+lineStrength, 0.01);

    // Don't have a real way to check result for individual wavelengths yet. I guess we could
    // try increasing the number of wavelengths and see if it converges to constant + line
    // profile value.
}

TEST_CASE("H line profile normalizations")
{
    // For each H line, make a line profile object and integrate it over a flat spectrum
    // here to check the normalization
    double rtol = .01;

    // Gather info about H lines
    HFromFiles hl;
    int numLines;
    Array lineFreqv;
    Array naturalLineWidthv;
    hl.lineInfo(numLines, lineFreqv, naturalLineWidthv);

    // Make flat spectrum
    double minFreq = Testing::defaultMinFreq;
    double maxFreq = Testing::defaultMaxFreq;
    Array someFreqs = {minFreq, (minFreq + maxFreq) / 2, maxFreq};
    // Array someFreqs = Testing::defaultCoarseFrequencyv();
    Spectrum flat(someFreqs, Array(1, someFreqs.size()));

    // Test the Lyman alpha line and another line
    for (int i : {0, 5})
    {
        double nu0 = lineFreqv[i];
        double halfWidth = naturalLineWidthv[i];

        // Try different gaussian broadenings
        for (double T : {0, 1, 10, 100, 1000, 10000})
        {
            double sigma_nu = nu0 * Functions::thermalVelocityWidth(T, Constant::HMASS) / Constant::LIGHT;
            LineProfile lp(nu0, sigma_nu, halfWidth);
            double norm = lp.integrateSpectrum(flat);
            bool line_norm_ok = TemplatedUtils::equalWithinTolerance(norm, 1., rtol);
            CHECK_MESSAGE(line_norm_ok, "norm of H line nr " << i << " is " << norm << " at " << T << " Kelvin");
        }
    }
}
