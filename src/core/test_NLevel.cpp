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
