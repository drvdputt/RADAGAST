#include "DoctestUtils.hpp"
#include "Ionization.hpp"
#include "Spectrum.hpp"
#include "SimpleHChemistry.hpp"
#include "SpeciesIndex.hpp"
#include "Testing.hpp"

TEST_CASE("SimpleHChemistry: compare exact solution of ionization")
{
	const double T = 10000;
	Array frequencyv = Testing::generateGeometricGridv(200, 1e11, 1e16);
	Array specificIntensityv = Testing::generateSpecificIntensityv(frequencyv, 25000, 10);
	Spectrum specificIntensity(frequencyv, specificIntensityv);

	SimpleHChemistry chemistry{};

	// no h2
	double kform = 0;
	double kdiss = 0;
	EVector kv = chemistry.rateCoeffv(T, specificIntensity, kdiss, kform);
	std::stringstream ss;
	ss << "Rate coeff: ionization, recombination, dissociation\n" << kv << '\n';
	std::string ratesMessage = ss.str();
	CAPTURE(ratesMessage);

	int ie = SpeciesIndex::index("e-");
	int ip = SpeciesIndex::index("H+");
	int iH = SpeciesIndex::index("H");
	int iH2 = SpeciesIndex::index("H2");

	EVector n0v(SpeciesIndex::size());
	n0v(ie) = 50;
	n0v(ip) = 50;
	n0v(iH) = 50;
	n0v(iH2) = 0;
	EVector nv = chemistry.solveBalance(kv, n0v);
	double f_network = nv(ip) / (nv(ip) + nv(iH));
	double f_exact = Ionization::solveBalance(nv(iH) + nv(ip), T, specificIntensity);
	CHECK_MESSAGE(TemplatedUtils::equalWithinTolerance(f_exact, f_network, 1e-6),
	              "exact ionization rate = "
	                              << f_exact << " while the one from chemical network is "
	                              << f_network << " (species vector is " << nv << ")");
}
