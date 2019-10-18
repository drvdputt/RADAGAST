#include "doctest.h"

#include "Constants.hpp"
#include "DoctestUtils.hpp"
#include "GasStruct.hpp"
#include "NLevel.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"
#include "TwoLevelHardcoded.hpp"

TEST_CASE("Test LevelCoefficients implementation two level subclass")
{
	TwoLevelHardcoded twolv{};
	EVector ev = twolv.ev();
	EVector gv = twolv.gv();
	EMatrix avv = twolv.avv();
	GasStruct gas;
	Spectrum specificIntensity;

	CHECK(twolv.numLv() == 2);

	int numLines;
	Array lineFreqv, naturalLineWidthv;
	twolv.lineInfo(numLines, lineFreqv, naturalLineWidthv);

	SUBCASE("lineInfo")
	{
		CHECK(numLines == 1);
		CHECK(numLines == twolv.numLv() - 1);
		CHECK(lineFreqv.size() == numLines);
		CHECK(naturalLineWidthv.size() == numLines);
		CHECK(lineFreqv[0] == (ev(1) - ev(0)) / Constant::PLANCK);
		CHECK(naturalLineWidthv[0] == avv(1, 0) / Constant::FPI);
	}

	SUBCASE("totalTransitionRatesvv")
	{
		EMatrix Tvv = twolv.totalTransitionRatesvv(specificIntensity, gas);
		CHECK(Tvv == avv);

		gas._speciesNv[SpeciesIndex::ine()] = 100;
		EMatrix cvv;
		Tvv = twolv.totalTransitionRatesvv(specificIntensity, gas, &cvv);
		CHECK(Tvv == avv + cvv);
		CHECK(cvv == twolv.cvv(gas));
	}
}
