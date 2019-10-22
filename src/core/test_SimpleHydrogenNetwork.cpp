#include "doctest.h"

#include "SimpleHydrogenNetwork.hpp"

#include "ChemistrySolver.hpp"
#include "Ionization.hpp"
#include "SpeciesIndex.hpp"
#include "Testing.hpp"

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
	n0v(ie) = 50;
	n0v(ip) = 50;
	n0v(iH) = 50;
	n0v(iH2) = 0;
	EVector nv = cs.solveBalance(kv, n0v);
	double f_network = nv(ip) / (nv(ip) + nv(iH));
	double f_exact = Ionization::solveBalance(nv(iH) + nv(ip), T, specificIntensity);
	CHECK_MESSAGE(TemplatedUtils::equalWithinTolerance(f_exact, f_network, 1e-6),
	              "exact ionization rate = "
	                              << f_exact << " while the one from chemical network is "
	                              << f_network << " (species vector is " << nv << ")");
}
