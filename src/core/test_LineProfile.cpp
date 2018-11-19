#include "doctest.h"

#include "Error.h"
#include "HydrogenFromFiles.h"
#include "LineProfile.h"
#include "NLevel.h"
#include "Testing.h"

namespace
{
LineProfile testLine()
{
	double c = 4.;
	double sg = 0.1;
	double hwl = 0.1;
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
		Error::fuzzyCheck("test value", spectrumv[i], manual, e);
	}
}

TEST_CASE("Test line profile integration over spectrum")
{
	auto lp = testLine();

	// Generate a flat spectrum of value base
	double base = 4.;
	auto s = testSpectrum(base);

	double integral = lp.integrateSpectrum(s);

	// This integral should be about equal to base, if gridpoints are chosen well
	double e = 1.e-2;
	Error::fuzzyCheck("test value", integral, base, e);
}

TEST_CASE("H line profile normalizations")
{
	// For each H line, make a line profile object and integrate it over a flat spectrum
	// here to check the normalization
	double rtol = .1;

	// Gather info about H lines
	NLevel hl(std::make_shared<HydrogenFromFiles>(), Constant::HMASS_CGS);
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

	for (int i = 3; i < numLines; i++)
	{
		double nu0 = lineFreqv[i];
		double halfWidth = naturalLineWidthv[i];

		// Choose a gaussian broadening of 100 K (0 causes nan)
		double thermalVelocity = sqrt(Constant::BOLTZMAN * 10000 / Constant::HMASS_CGS);
		double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;

		LineProfile lp(nu0, sigma_nu, halfWidth);
		double norm = lp.integrateSpectrum(flat);
		std::cout << "H line nr " << i << std::endl;
		bool line_norm_ok = TemplatedUtils::equalWithinTolerance(norm, 1., rtol);
		WARN_MESSAGE(line_norm_ok, "norm of H line nr " << i << " is " << norm);
	}
}
