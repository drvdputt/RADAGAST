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

	EVector n0v(4);
	n0v(ie) = 0;
	n0v(ip) = 0;
	n0v(iH) = 100;
	n0v(iH2) = 0;

	EVector nv = cs.solveBalance(kv, n0v);
	double ionizedFraction =
	                Ionization::solveBalance(nv(iH) + nv(ip), T, specificIntensity);

	std::cout << "Compare with ionized fraction calculation: " << std::endl;
	std::cout << "f = " << ionizedFraction << std::endl;
	CHECK(TemplatedUtils::equalWithinTolerance(nv(ip) / (nv(ip) + nv(iH)), ionizedFraction, 0.00));
}

// TEST_CASE("single species which dissolves into nothing")
// {
// }
