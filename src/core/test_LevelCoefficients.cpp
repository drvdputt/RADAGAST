#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DoctestUtils.hpp"
#include "LevelCoefficients.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"
#include "TwoLevelHardcoded.hpp"
#include "doctest.h"

using namespace GasModule;

TEST_CASE("Test LevelCoefficients implementation two level subclass")
{
    TwoLevelHardcoded twolv{};
    EVector ev = twolv.ev();
    EVector gv = twolv.gv();
    EMatrix avv = twolv.avv();

    SpeciesIndex spindex({"e-"});
    SpeciesVector sv(&spindex);
    CollisionParameters gas(500, sv);

    Spectrum meanIntensity;

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
}
