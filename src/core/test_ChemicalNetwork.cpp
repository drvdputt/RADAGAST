#include "doctest.h"

#include "ChemistrySolver.h"
#include "IonizationBalance.h"
#include "SimpleHydrogenNetwork.h"
#include "SpeciesIndex.h"
#include "Testing.h"

TEST_CASE("testing chemistry solver by comparing with exact solution of ionization")
{
	const double T = 10000;
	Array frequencyv = Testing::generateGeometricGridv(200, 1e11, 1e16);
	Array specificIntensityv = Testing::generateSpecificIntensityv(frequencyv, 25000, 10);
	Spectrum specificIntensity(frequencyv, specificIntensityv);

	ChemistrySolver cs(std::make_unique<SimpleHydrogenNetwork>());

	// Formation and dissociation rates should come from somewhere else
	double kform = 0;
	double kdiss = 0;

	EVector kv = cs.chemicalNetwork()->rateCoeffv(T, specificIntensity, kdiss, kform);
	std::cout << "Rate coeff: ionization, recombination, dissociation" << std::endl
	          << kv << std::endl;

	int ie = SpeciesIndex::index("e-");
	int ip = SpeciesIndex::index("H+");
	int iH = SpeciesIndex::index("H");
	int iH2 = SpeciesIndex::index("H2");

	EVector n0v(SpeciesIndex::size());
	n0v(ie) = 0;
	n0v(ip) = 0;
	n0v(iH) = 100;
	n0v(iH2) = 0;

	EVector nv = cs.solveBalance(kv, n0v);
	double ionizedFraction =
	                Ionization::solveBalance(nv(iH) + nv(ip), T, specificIntensity);

	std::cout << "Compare with ionized fraction calculation: " << std::endl;
	std::cout << "f = " << ionizedFraction << std::endl;
	CHECK(TemplatedUtils::equalWithinTolerance(nv(ip) / (nv(ip) + nv(iH)), ionizedFraction,
	                                           0.00));
}

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
		double c = 10.;
		double d = 15.;
		CHECK(solve(c, d) == c / d);
	}

	SUBCASE("only destruction, should tend to zero") { CHECK(solve(0, 1.) == 0.); }

	SUBCASE("only creation, should go to infinity, technically") {
		double solution = solve(0., 1.);
		WARN_MESSAGE(std::isinf(solve(0, 1.)), "solution is " << solution);
	}
}

TEST_CASE("Combine and dissociate")
{
	// two species (e.g. H and H2), two reactions
	EMatrix r(2,2);
	EMatrix p(2,2);

	// combination
	r.col(0) << 2, 0; // 2 H consumed
	p.col(0) << 0, 1; // 1 H2 produced

	// dissociation
	r.col(1) << 0, 1; // 1 H2 consumed
	p.col(1) << 2, 0; // 2 H2 produced

	// conservation equation
	// 1 equation, 2 species
	EMatrix c(1, 2);
	c << 1, 2; // 1 * nH + 2 * nH2 is conserved (number of protons)

	ChemistrySolver cs(r, p, c);

	double nH0 = 10;
	double nH20 = 20;
	double Ntotal = nH0 + 2 * nH20;
	EVector n0v(2);
	n0v << nH0, nH20;

	double eps = 1.e-15;
	SUBCASE("balance")
	{
		double kform = 1.;
		double kdiss = 0.2;
		EVector kv(2);
		kv << kform, kdiss;

		// Exact solution
		double nH_exact = (kdiss / 2 - std::sqrt(kdiss * kdiss / 4 + 2 * kdiss * kform * Ntotal)) / (-2 * kform);
		double nH2_exact = (Ntotal - nH_exact) / 2;

		// General algorithm
		EVector nv = cs.solveBalance(kv, n0v);
		CHECK(nv(0) <= nH_exact + eps);
		CHECK(nv(0) >= nH_exact - eps);
		CHECK(nv(1) == nH2_exact);
	}
	SUBCASE("only formation")
	{
		double kform = 1.;
		double kdiss = 0;
		EVector kv(2);
		kv << kform, kdiss;

		EVector nv = cs.solveBalance(kv, n0v);
		// All H should disappear, and be transformed into H2
		CHECK(nv(0) >= 0.);
		CHECK(nv(0) < eps);
		CHECK(nv(1) == Ntotal / 2.);
	}
	SUBCASE("only dissociation")
	{
		double kform = 0;
		double kdiss = 1.;
		EVector kv(2);
		kv << kform, kdiss;

		EVector nv = cs.solveBalance(kv, n0v);
		// All H2 should disappear, and be transformed into H
		CHECK(nv(0) == Ntotal);
		CHECK(nv(1) >= 0.);
		CHECK(nv(1) < eps);
	}
}
