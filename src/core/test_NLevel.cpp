#include "doctest.h"

#include "Constants.hpp"
#include "DoctestUtils.hpp"
#include "GasStruct.hpp"
#include "NLevel.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"
#include "TwoLevelHardcoded.hpp"

TEST_CASE("Test NLevel implementation using two level system")
{
	TwoLevelHardcoded twolv{};
	EVector ev = twolv.ev();
	EVector gv = twolv.gv();
	EMatrix avv = twolv.avv();
	GasStruct gas;
	Spectrum specificIntensity;

	NLevel nlv{std::make_shared<TwoLevelHardcoded>(twolv), Constant::HMASS_CGS * 12};
	CHECK(nlv.numLv() == 2);

	int numLines;
	Array lineFreqv, naturalLineWidthv;
	nlv.lineInfo(numLines, lineFreqv, naturalLineWidthv);

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

		double nu = lineFreqv[0];
		double c = Constant::LIGHT;
		double h = Constant::PLANCK;
		double expectedRatio = 2 * h * nu * nu * nu * s.nv(1) / c / c / (s.nv(0) * gv(1) / gv(0) - s.nv(1));

		// Check integrated and individual emissivity / opacity ratios, to check for
		// conservation when integrating.
		// double integratedEm = TemplatedUtils::integrate<double>(eFrequencyv, emissivityv);
		// double integratedOp = TemplatedUtils::integrate<double>(eFrequencyv, opacityv);
		// CHECK(integratedEm / integratedOp == expectedRatio);
		double eps = 1e-15;
		for (size_t iFreq = 0; iFreq < eFrequencyv.size(); iFreq++)
		{
			double op = opacityv[iFreq];
			double em = emissivityv[iFreq];
			if (op != 0)
				DoctestUtils::checkTolerance("individual em / op ratio", em / op, expectedRatio, eps);
		}

		// heating - cooling should be zero in LTE
		double heat = nlv.heating(s);
		double cool = nlv.cooling(s);
		CHECK(heat - cool == 0);
	}
}
