#include "BigH2Model.hpp"
#include "CollisionParameters.hpp"
#include "DoctestUtils.hpp"
#include "H2Data.hpp"
#include "SpecialFunctions.hpp"
#include "SpeciesIndex.hpp"
#include "Testing.hpp"

TEST_CASE("H2-specific algorithm")
{
    H2Data hff(99, 99);
    BigH2Model bigh2(&hff);

    EVector ev = hff.ev();

    double T = 500;

    SpeciesIndex spindex(SpeciesIndex::e_p_H_H2);
    EVector zerov = EVector::Zero(hff.numLv());
    Array frequencyv = Testing::generateGeometricGridv();

    double epsFrac = 1e-2;

    SpeciesVector sv(&spindex);
    auto makeCP = [&](double ne, double np, double nH, double nH2) {
        sv.setDensities(spindex.linearCombination(SpeciesIndex::e_p_H_H2, {ne, np, nH, nH2}));
        return CollisionParameters(T, sv);
    };

    SUBCASE("no radiation, high density, should go to LTE")
    {
        double n = 1e8;
        const CollisionParameters cp = makeCP(100, 100, 100, n);

        Array specificIntensityv(frequencyv.size());
        Spectrum specificIntensity(frequencyv, specificIntensityv);

        bigh2.solve(n, cp, specificIntensity);
        EVector nv_lte = n * hff.solveBoltzmanEquations(T);

        const EVector& nv_sol = bigh2.levelSolution()->nv();

        // This test seems to work reasonably for the first three levels
        for (size_t i = 0; i < std::min<int>(16, nv_lte.size()); i++)
            DoctestUtils::checkTolerance("level " + std::to_string(i) + " from solver vs LTE", nv_sol(i), nv_lte(i),
                                         epsFrac);

        // Check if the total density is correct
        double eps = 1e-15;
        DoctestUtils::checkTolerance("sum of s0.nv", nv_sol.sum(), n, eps);
        DoctestUtils::checkTolerance("sum of sLTE.nv", nv_lte.sum(), n, eps);
    }

    SUBCASE("no radiation, low density, should go to ground state")
    {
        double n = 1e-15;
        double eps = 1e-2;

        const CollisionParameters cp = makeCP(0, 0, 0, n);

        Array specificIntensityv(frequencyv.size());
        Spectrum specificIntensity(frequencyv, specificIntensityv);
        bigh2.solve(n, cp, specificIntensity);
        const EVector& nv = bigh2.levelSolution()->nv();

        // Check if some individual levels are indeed close to 0
        for (size_t i = 2; i < std::min<int>(16, nv.size()); i++)
            DoctestUtils::checkRange("level " + std::to_string(i), nv(i), 0, n * eps);

        // Check if almost all density is in lowest (two) level(s) (ortho+para)
        DoctestUtils::checkTolerance("ground level(s)", nv(0) + nv(1), n, eps);

        CAPTURE("ground level (energy " << ev(0) << ") " << nv(0) << " \n level 1 (energy " << ev(1) << ") " << nv(1)
                                        << " \n level 3 (energy " << ev(3) << ") " << nv(3) << "\n");

        // Check if all but the lowest two level populations are small compared to total
        CHECK((nv.tail(nv.size() - 2).array() < eps * n).all());
    }

    // I need a way to check that the effect of radiation on level populations is handles
    // correctly. The following only seems to work when the initial guess is LTE (which is
    // standard now).
    SUBCASE("Blackbody radiation, no collisions, should also go to LTE?")
    {
        double n = 1;
        const CollisionParameters cp = makeCP(0, 0, 0, n);
        Array specificIntensityv(frequencyv.size());
        for (size_t i = 0; i < frequencyv.size(); i++)
            specificIntensityv[i] = SpecialFunctions::planck(frequencyv[i], T);
        Spectrum specificIntensity(frequencyv, specificIntensityv);
        bigh2.solve(n, cp, specificIntensity);
        EVector nv_lte = n * hff.solveBoltzmanEquations(T);
        const EVector& nv_sol = bigh2.levelSolution()->nv();

        // This test seems to work reasonably for the first bunch of levels (fails at
        // 16, and produces warning)
        for (size_t i = 0; i < std::min<int>(16, nv_lte.size()); i++)
            DoctestUtils::checkTolerance("level " + std::to_string(i) + " from solver vs LTE", nv_sol(i), nv_lte(i),
                                         epsFrac, true);

        // Check if the total density is correct
        double eps = 1e-15;
        DoctestUtils::checkTolerance("sum of s0.nv", nv_sol.sum(), n, eps);
        DoctestUtils::checkTolerance("sum of sLTE.nv", nv_lte.sum(), n, eps);
    }
}
