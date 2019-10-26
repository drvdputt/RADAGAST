#include "Chemistry.hpp"
#include "DoctestUtils.hpp"
#include "SpeciesIndex.hpp"

// Note that the species in this test are purely fictional. I'm just using existing
// names/indices from the SpeciesIndex until the latter is reworked.

TEST_CASE("single species with creation and destruction")
{
	// one species, two reactions
	Chemistry chemistry;
	chemistry.addReaction("creation", {}, {}, {"e-"}, {1});
	chemistry.addReaction("destruction", {"e-"}, {1}, {}, {});
	chemistry.prepareCoefficients();

	auto solve = [&](double creationRate, double destructionCoeff) -> double {
		EVector k(2);
		k(chemistry.reactionIndex("creation")) = creationRate;
		k(chemistry.reactionIndex("destruction")) = destructionCoeff;
		k << creationRate, destructionCoeff;

		EVector n0 = EVector::Zero(SpeciesIndex::size());
		n0(SpeciesIndex::ine()) = 100.;
		EVector n = chemistry.solveBalance(k, n0);
		return n(SpeciesIndex::ine());
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
	DoctestUtils::checkRange("single species, only destruction, should tend to zero",
	                         onlyDestruction, 0, eps);
}

TEST_CASE("Combine and dissociate")
{
	// two species (e.g. H and H2), two reactions
	Chemistry chemistry;
	chemistry.addReaction("combine", {"H"}, {2}, {"H2"}, {1});
	chemistry.addReaction("dissociate", {"H2"}, {1}, {"H"}, {2});
	chemistry.prepareCoefficients();

	// Try different initial conditions (TEST_CASE is restarted for every subcase)
	double nH0;
	double nH20;
	SUBCASE("start with only H")
	{
		nH0 = 30;
		nH20 = 0;
	}
	SUBCASE("start with only H2")
	{
		nH0 = 0;
		nH20 = 15;
	}
	SUBCASE("start with both H and H2")
	{
		nH0 = 10;
		nH20 = 10;
	}
	SUBCASE("small large")
	{
		nH0 = 1e-10;
		nH20 = 1000;
	}
	SUBCASE("large small")
	{
		nH0 = 1e6;
		nH20 = 1e-6;
	}
	CAPTURE(nH0);
	CAPTURE(nH20);

	// For each of these subcases, we do the following tests
	double Ntotal = nH0 + 2 * nH20;
	EVector n0v = EVector::Zero(SpeciesIndex::size());
	n0v(SpeciesIndex::inH()) = nH0;
	n0v(SpeciesIndex::inH2()) = nH20;

	{ // Both formation and destruction
		double eps = 1.e-13;
		CAPTURE("Balance for " << nH0 << " " << nH20);
		double kform = 1.;
		double kdiss = 0.2;
		EVector kv(2);
		kv(chemistry.reactionIndex("combine")) = kform;
		kv(chemistry.reactionIndex("dissociate")) = kdiss;

		// Exact solution
		double nH_exact = (kdiss / 2 -
		                   std::sqrt(kdiss * kdiss / 4 + 2 * kdiss * kform * Ntotal)) /
		                  (-2 * kform);
		double nH2_exact = (Ntotal - nH_exact) / 2;

		// General algorithm
		EVector nv = chemistry.solveBalance(kv, n0v);
		DoctestUtils::checkTolerance("nH (should be analytic solution)",
		                             nv(SpeciesIndex::inH()), nH_exact, eps);
		DoctestUtils::checkTolerance("nH2 and nH2_exact", nv(SpeciesIndex::inH2()),
		                             nH2_exact, eps);
	}

	{ // Only formation
		double eps = 1.e-10;
		double kform = 1.;
		double kdiss = 0;
		EVector kv(2);
		kv(chemistry.reactionIndex("combine")) = kform;
		kv(chemistry.reactionIndex("dissociate")) = kdiss;

		EVector nv = chemistry.solveBalance(kv, n0v);
		// All H should disappear, and be transformed into H2
		DoctestUtils::checkRange("nH (should disappear)", nv(SpeciesIndex::inH()), 0.,
		                         eps);
		DoctestUtils::checkTolerance("nH2", nv(SpeciesIndex::inH2()), Ntotal / 2., eps);
	}

	{ // Only dissociation
		double eps = 1.e-10;
		double kform = 0;
		double kdiss = 1.;
		EVector kv(2);
		kv(chemistry.reactionIndex("combine")) = kform;
		kv(chemistry.reactionIndex("dissociate")) = kdiss;

		EVector nv = chemistry.solveBalance(kv, n0v);
		// All H2 should disappear, and be transformed into H
		DoctestUtils::checkRange("nH2 (should disappear)", nv(SpeciesIndex::inH2()), 0.,
		                         eps);
		DoctestUtils::checkTolerance("nH2 (should equal total)",
		                             nv(SpeciesIndex::inH()), Ntotal, eps);
	}
}
