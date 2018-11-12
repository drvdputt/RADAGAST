#include "doctest.h"

#include "Constants.h"
#include "GasStruct.h"
#include "NLevel.h"
#include "SpeciesIndex.h"
#include "Spectrum.h"
#include "Testing.h"
#include "TwoLevelHardcoded.h"

TEST_CASE("Test NLevel implementation using two level system")
{
	TwoLevelHardcoded twolv{};
	EVector ev = twolv.ev();
	EMatrix avv = twolv.avv();
	GasStruct gas;
	Spectrum specificIntensity;

	NLevel nlv{std::make_shared<TwoLevelHardcoded>(twolv), Constant::HMASS_CGS * 12};

	SUBCASE("lineInfo")
	{
		int numLines;
		Array lineFreqv, naturalLineWidthv;
		nlv.lineInfo(numLines, lineFreqv, naturalLineWidthv);
		CHECK(numLines == 1);
		CHECK(numLines == twolv.numLv() - 1);
		CHECK(lineFreqv.size() == numLines);
		CHECK(naturalLineWidthv.size() == numLines);
		CHECK(lineFreqv[0] == (ev(1) - ev(0)) / Constant::PLANCK);
		CHECK(naturalLineWidthv[0] == avv(1, 0) / Constant::FPI);
	}

	SUBCASE("totalTransitionRatesvv")
	{
		EMatrix Tvv = nlv.totalTransitionRatesvv(specificIntensity, gas);
		CHECK(Tvv == avv);

		EMatrix cvv;
		Tvv = nlv.totalTransitionRatesvv(specificIntensity, gas, &cvv);
		CHECK(Tvv == avv + cvv);
		CHECK(cvv == twolv.cvv(gas));
	}

	SUBCASE("LTE")
	{
		NLevel::Solution s = nlv.solveLTE(1, gas);
		Array eFrequencyv = Testing::generateGeometricGridv();
		Array emissivityv = nlv.emissivityv(s, eFrequencyv);
		Array opacityv = nlv.opacityv(s, eFrequencyv);
	}
}

// #define REPORT_LINE_QUALITY
#ifdef REPORT_LINE_QUALITY
// Full integral of the line profile, to check the discretization of the output
// (emission) grid.
double maxNorm = 0, minNorm = 1e9;
forActiveLinesDo([&](size_t upper, size_t lower) {
	auto lp = lineProfile(upper, lower, s.T, s.cvv);

	Array lpv(specificIntensity.frequencyv().size());
	for (size_t i = 0; i < lpv.size(); i++)
		lpv[i] = lp(specificIntensity.frequencyv()[i]);

	double norm = TemplatedUtils::integrate<double>(specificIntensity.frequencyv(), lpv);
	DEBUG("on the input intensity grid, line  " << upper << " --> " << lower << " has norm "
	                                            << norm << endl);
	maxNorm = max(norm, maxNorm);
	minNorm = min(norm, minNorm);
});
DEBUG("Max profile norm = " << maxNorm << endl);
DEBUG("Min profile norm = " << minNorm << endl);
maxNorm = 0;
minNorm = 1e9; // any large number will do
// Integral over contant spectrum, using the LineProfile class. An integration
// grid is chosen internally, and we check the quality of it here (at least for
// now.)
double minFreq = specificIntensity.freqMin();
double maxFreq = specificIntensity.freqMax();
Array someFreqs = {minFreq, (minFreq + maxFreq) / 2, maxFreq};
Spectrum flat(someFreqs, Array(1, someFreqs.size()));
forActiveLinesDo([&](size_t upper, size_t lower) {
	auto lp = lineProfile(upper, lower, s);
	double norm = lp.integrateSpectrum(flat);
	DEBUG("LineProfile " << upper << " --> " << lower << " has norm " << norm << endl);
	maxNorm = max(norm, maxNorm);
	minNorm = min(norm, minNorm);
});
DEBUG("Max profile norm = " << maxNorm << endl);
DEBUG("Min profile norm = " << minNorm << endl);
#endif /* REPORT_LINE_QUALITY */
