#include "CollisionParameters.hpp"
#include "DoctestUtils.hpp"
#include "LevelSolution.hpp"
#include "Testing.hpp"
#include "TwoLevelHardcoded.hpp"
#include "doctest.h"

using namespace GasModule;

TEST_CASE("Test LevelSolver using two-level LTE ")
{
    double T = 500;
    double n = 1000;
    double ne = n;
    SpeciesIndex spindex({"H", "e-"});
    SpeciesVector sv(&spindex);
    sv.setNe(ne);
    CollisionParameters cp(T, sv);

    TwoLevelHardcoded twolv;
    LevelSolution s(&twolv);
    s.updateRates(cp);
    s.setNv(n * twolv.solveBoltzmanEquations(T));

    // Deliberately use very low resolution
    Array eFrequencyv = Testing::generateGeometricGridv(5, Testing::defaultMinFreq, Testing::defaultMaxFreq);

    Array emissivityv = s.emissivityv(eFrequencyv);
    Array opacityv = s.opacityv(eFrequencyv);

    int numLines;
    Array lineFreqv, naturalLineWidthv;
    twolv.lineInfo(numLines, lineFreqv, naturalLineWidthv);
    double nu = lineFreqv[0];

    double c = Constant::LIGHT;
    double h = Constant::PLANCK;
    double expectedRatio =
        2 * h * nu * nu * nu * s.nv()(1) / c / c / (s.nv()(0) * twolv.gv()(1) / twolv.gv()(0) - s.nv()(1));

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
        if (op != 0) DoctestUtils::checkTolerance("individual em / op ratio", em / op, expectedRatio, eps);
    }

    // heating - cooling should be zero in LTE
    double heat = s.netHeating();
    CAPTURE(heat);
    CHECK(abs(heat) < 1.e-30);
}
