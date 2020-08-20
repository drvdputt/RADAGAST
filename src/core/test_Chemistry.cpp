#include "Chemistry.hpp"
#include "DoctestUtils.hpp"
#include "SpeciesIndex.hpp"

using namespace RADAGAST;

TEST_CASE("single species with creation and destruction")
{
    // one species, two reactions
    Chemistry chemistry;
    chemistry.registerSpecies({"e-"});
    chemistry.addReaction("creation", {}, {}, {"e-"}, {1});
    chemistry.addReaction("destruction", {"e-"}, {1}, {}, {});
    chemistry.prepareCoefficients();

    auto solve = [&](double creationRate, double destructionCoeff) -> double {
        EVector k(2);
        k(chemistry.reactionIndex("creation")) = creationRate;
        k(chemistry.reactionIndex("destruction")) = destructionCoeff;
        k << creationRate, destructionCoeff;

        SpeciesVector sv0(&chemistry.speciesIndex());
        sv0.setNe(100.);
        SpeciesVector sv(&chemistry.speciesIndex());
        sv.setDensities(chemistry.solveBalance(k, sv0.speciesNv()));
        return sv.ne();
    };

    SUBCASE("two thirds ratio")
    {
        // equilibrium means:
        // creation - N * destruction = 0 --> N = creation / destruction
        double cr = 10.;
        double ds = 15.;
        CHECK(solve(cr, ds) == cr / ds);
    }

    double eps = 1e-15;
    double onlyDestruction = solve(0, 1.);
    DoctestUtils::checkRange("single species, only destruction, should tend to zero", onlyDestruction, 0, eps);
}

TEST_CASE("Combine and dissociate")
{
    // two species (e.g. H and H2), two reactions. Also tests abstraction by calling them A
    // and B.
    Chemistry chemistry;
    chemistry.registerSpecies({"A", "B"});
    chemistry.addReaction("combine", {"A"}, {2}, {"B"}, {1});
    chemistry.addReaction("dissociate", {"B"}, {1}, {"A"}, {2});
    chemistry.prepareCoefficients();

    // Try different initial conditions (TEST_CASE is restarted for every subcase)
    double nA;
    double nB;
    SUBCASE("start with only A")
    {
        nA = 30;
        nB = 0;
    }
    SUBCASE("start with only B")
    {
        nA = 0;
        nB = 15;
    }
    SUBCASE("start with both A and B")
    {
        nA = 10;
        nB = 10;
    }
    SUBCASE("small large")
    {
        nA = 1e-10;
        nB = 1000;
    }
    SUBCASE("large small")
    {
        nA = 1e6;
        nB = 1e-6;
    }
    CAPTURE(nA);
    CAPTURE(nB);

    // For each of these subcases, we do the following tests
    double Ntotal = nA + 2 * nB;
    EVector n0v = EVector::Zero(chemistry.numSpecies());
    n0v(chemistry.speciesIndex().index("A")) = nA;
    n0v(chemistry.speciesIndex().index("B")) = nB;

    SpeciesVector sv(&chemistry.speciesIndex());

    {  // Both formation and destruction
        double eps = 1.e-13;
        CAPTURE("Balance for " << nA << " " << nB);
        double kform = 1.;
        double kdiss = 0.2;
        EVector kv(2);
        kv(chemistry.reactionIndex("combine")) = kform;
        kv(chemistry.reactionIndex("dissociate")) = kdiss;

        // Exact solution
        double nA_exact = (kdiss / 2 - std::sqrt(kdiss * kdiss / 4 + 2 * kdiss * kform * Ntotal)) / (-2 * kform);
        double nB_exact = (Ntotal - nA_exact) / 2;

        // General algorithm
        sv.setDensities(chemistry.solveBalance(kv, n0v));
        DoctestUtils::checkTolerance("nA (should be analytic solution)", sv.nSpecies("A"), nA_exact, eps);
        DoctestUtils::checkTolerance("nB and nB_exact", sv.nSpecies("B"), nB_exact, eps);
    }

    {  // Only formation
        double eps = 1.e-10;
        double kform = 1.;
        double kdiss = 0;
        EVector kv(2);
        kv(chemistry.reactionIndex("combine")) = kform;
        kv(chemistry.reactionIndex("dissociate")) = kdiss;

        sv.setDensities(chemistry.solveBalance(kv, n0v));
        // All H should disappear, and be transformed into H2
        DoctestUtils::checkRange("nA (should disappear)", sv.nSpecies("A"), 0., eps);
        DoctestUtils::checkTolerance("nB", sv.nSpecies("B"), Ntotal / 2., eps);
    }

    {  // Only dissociation
        double eps = 1.e-10;
        double kform = 0;
        double kdiss = 1.;
        EVector kv(2);
        kv(chemistry.reactionIndex("combine")) = kform;
        kv(chemistry.reactionIndex("dissociate")) = kdiss;

        sv.setDensities(chemistry.solveBalance(kv, n0v));
        // All H2 should disappear, and be transformed into H
        DoctestUtils::checkRange("nH2 (should disappear)", sv.nSpecies("B"), 0., eps);
        DoctestUtils::checkTolerance("nH2 (should equal total)", sv.nSpecies("A"), Ntotal, eps);
    }
}
