#include "doctest.h"

#include "DoctestUtils.hpp"
#include "HFromFiles.hpp"
#include "LineProfile.hpp"
#include "LevelCoefficients.hpp"
#include "Testing.hpp"

namespace
{
LineProfile testLine()
{
	double c = 4.;
	double sg = 0.1;
	double hwl = 0.1;
	return LineProfile(c, sg, hwl);
}

LineProfile mostlyGaussianLine()
{
	double c = 4.;
	double sg = 0.2;
	double hwl = 0.01;
	return LineProfile(c, sg, hwl);
}

Spectrum testSpectrum(double base)
{
	size_t numFreq = 100;
	Array frequencyv(numFreq);
	double fmin = .01;
	double fmax = 10.;
	double step = (fmax - fmin) / numFreq;
	for (size_t i = 0; i < numFreq; i++)
		frequencyv[i] = fmin + i * step;
	Array spectrumv(base, numFreq);
	return Spectrum(frequencyv, spectrumv);
}
} // namespace

TEST_CASE("Test line profile addition to spectrum")
{
	auto lp = testLine();

	double base = 1.;
	auto s = testSpectrum(base);
	Array frequencyv = s.frequencyv();
	Array spectrumv = s.valuev();

	double factor = 3.;
	lp.addToSpectrum(frequencyv, spectrumv, factor);

	double e = 1.e-15;
	for (size_t i = 0; i < frequencyv.size(); i++)
	{
		double linevalue = factor * lp(frequencyv[i]);
		double manual = base + linevalue;
		DoctestUtils::checkTolerance("test value", spectrumv[i], manual, e);
	}
}

TEST_CASE("Test line profile integration over spectrum")
{
	// Generate a flat spectrum of value base
	double base = 4.;
	auto s = testSpectrum(base);

	SUBCASE("mixed guassian/lorentzian line")
	{
		auto lp = testLine();
		double integral = lp.integrateSpectrum(s, 0, "guasslorentz-points.dat");
		// This integral should be about equal to base, if gridpoints are chosen well
		double e = 1.e-2;
		DoctestUtils::checkTolerance("test value", integral, base, e, true);
	}

	SUBCASE("mostly gaussian line")
	{
		auto lp = mostlyGaussianLine();
		double integral = lp.integrateSpectrum(s, 0, "gauss-points.dat");
		// This integral should be about equal to base, if gridpoints are chosen well
		double e = 1.e-2;
		DoctestUtils::checkTolerance("test value", integral, base, e);
	}
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
		for (double T : {1, 10, 100, 1000, 10000})
		{
			double thermalVelocity =
			                sqrt(Constant::BOLTZMAN * T / Constant::HMASS_CGS);
			double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;

			LineProfile lp(nu0, sigma_nu, halfWidth);
			double norm = lp.integrateSpectrum(flat);
			bool line_norm_ok =
			                TemplatedUtils::equalWithinTolerance(norm, 1., rtol);
			CHECK_MESSAGE(line_norm_ok, "norm of H line nr " << i << " is " << norm
			                                                << " at " << T
			                                                << " Kelvin");
		}
	}
}
