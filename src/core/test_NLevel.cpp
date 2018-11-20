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
	CHECK(nlv.numLv() == 2);

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

	SUBCASE("Test functions that need solution using LTE ")
	{
		NLevel::Solution s = nlv.solveLTE(1, gas);

		// Deliberately use very low resolution
		Array eFrequencyv = Testing::generateGeometricGridv(5, Testing::defaultMinFreq,
		                                                    Testing::defaultMaxFreq);

		Array emissivityv = nlv.emissivityv(s, eFrequencyv);
		Array opacityv = nlv.opacityv(s, eFrequencyv);
		Array lineEmv = nlv.lineEmissivityv(s, eFrequencyv);
		Array lineOpv = nlv.lineOpacityv(s, eFrequencyv);

		// Since we are working with the base class, there should be only line emission
		for (size_t iFreq = 0; iFreq < eFrequencyv.size(); iFreq++)
		{
			CHECK(emissivityv[iFreq] == lineEmv[iFreq]);
			CHECK(opacityv[iFreq] == lineOpv[iFreq]);
		}

		// For each frequency, the ratio should be equal to the emissivity/opacity
		// factor of the line
		Array em_op = lineEmv / lineOpv;
		WARN_MESSAGE(em_op[0] == 1,
		             "Check ratio of emissivity/opacity here, not implemented yet");

		// heating - cooling should be zero in LTE
		double heat = nlv.heating(s);
		double cool = nlv.cooling(s);
		CHECK(heat - cool == 0);
	}
}
