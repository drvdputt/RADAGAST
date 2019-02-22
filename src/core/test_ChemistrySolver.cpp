#include "ChemistrySolver.h"
#include "DoctestUtils.h"

TEST_CASE("single species with creation and destruction")
{
	// one species, two reactions
	EMatrix r(1, 2);
	EMatrix p(1, 2);

	// Creation (element 0, reaction 0)
	r(0, 0) = 0; // zero consumed
	p(0, 0) = 1; // one created

	// Destruction (element 0, reaction 1)
	r(0, 1) = 1; // one consumed
	p(0, 1) = 0; // nothing created

	// No conservation equations (good edge case to test)
	EMatrix c(0, 1);

	// Use this constructor to create a solver for this abstract chemical network (reaction
	// rates will have to be provided manually)
	ChemistrySolver cs(r, p, c);

	auto solve = [&](double creationRate, double destructionCoeff) -> double {
		EVector k(2);
		k << creationRate, destructionCoeff;

		EVector n0(1);
		n0(0) = 100.;
		EVector n = cs.solveBalance(k, n0);
		return n(0);
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
	EMatrix r(2, 2);
	EMatrix p(2, 2);

	// combination
	r.col(0) << 2, 0; // 2 H consumed
	p.col(0) << 0, 1; // 1 H2 produced

	// dissociation
	r.col(1) << 0, 1; // 1 H2 consumed
	p.col(1) << 2, 0; // 2 H produced

	// conservation equation
	// 1 equation, 2 species
	EMatrix c(1, 2);
	c << 1, 2; // 1 * nH + 2 * nH2 is conserved (number of protons)

	ChemistrySolver cs(r, p, c);

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
	EVector n0v(2);
	n0v << nH0, nH20;

	{ // Both formation and destruction
		double eps = 1.e-15;
		std::cout << "balance for " << nH0 << " " << nH20 << std::endl;
		double kform = 1.;
		double kdiss = 0.2;
		EVector kv(2);
		kv << kform, kdiss;

		// Exact solution
		double nH_exact = (kdiss / 2 -
		                   std::sqrt(kdiss * kdiss / 4 + 2 * kdiss * kform * Ntotal)) /
		                  (-2 * kform);
		double nH2_exact = (Ntotal - nH_exact) / 2;

		// General algorithm
		EVector nv = cs.solveBalance(kv, n0v);
		DoctestUtils::checkTolerance("nH (should be analytic solution)", nv(0),
		                             nH_exact, eps);
		DoctestUtils::checkTolerance("nH2 and nH2_exact", nv(1), nH2_exact, eps);
	}

	{ // Only formation
		double eps = 1.e-10;
		double kform = 1.;
		double kdiss = 0;
		EVector kv(2);
		kv << kform, kdiss;

		EVector nv = cs.solveBalance(kv, n0v);
		// All H should disappear, and be transformed into H2
		DoctestUtils::checkRange("nH (should disappear)", nv(0), 0., eps);
		DoctestUtils::checkTolerance("nH2", nv(1), Ntotal / 2., eps);
	}

	{ // Only dissociation
		double eps = 1.e-10;
		double kform = 0;
		double kdiss = 1.;
		EVector kv(2);
		kv << kform, kdiss;

		EVector nv = cs.solveBalance(kv, n0v);
		// All H2 should disappear, and be transformed into H
		DoctestUtils::checkRange("nH2 (should disappear)", nv(1), 0., eps);
		DoctestUtils::checkTolerance("nH2 (should equal total)", nv(0), Ntotal, eps);
	}
}
