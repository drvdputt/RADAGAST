#include "doctest.h"

#include "GasStruct.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "SpeciesIndex.h"
#include "Testing.h"

TEST_CASE("H2-specific algorithm")
{
	H2FromFiles hff(2, 1);
	H2Levels h2l(std::make_shared<const H2FromFiles>(hff));

	double T = 500;

	EVector speciesNv = EVector::Zero(SpeciesIndex::size());
	EVector zerov = EVector::Zero(h2l.numLv());
	Array frequencyv = Testing::generateGeometricGridv();

	double epsFrac = 1e-4;

	SUBCASE("no radiation, high density, should go to LTE")
	{
		double n = 1e5;
		speciesNv(SpeciesIndex::ine()) = 100;
		speciesNv(SpeciesIndex::inp()) = 100;
		speciesNv(SpeciesIndex::inH()) = 100;
		speciesNv(SpeciesIndex::inH2()) = n;
		const GasStruct gas(T, speciesNv);

		Array specificIntensityv(frequencyv.size());
		Spectrum specificIntensity(frequencyv, specificIntensityv);
		NLevel::Solution s0 = h2l.solveBalance(n, specificIntensity, zerov, zerov, gas);
		NLevel::Solution sLTE = h2l.solveLTE(n, specificIntensity, gas);

		// This test seems to work reasonable for the first three levels
		for (size_t i = 0; i < std::min<int>(3, s0.nv.size()); i++)
		{
			CHECK_MESSAGE(TemplatedUtils::equalWithinTolerance(s0.nv(i), sLTE.nv(i),
			                                                   epsFrac),
			              "level " << i << " from solver: " << s0.nv(i)
			                       << ", from LTE:" << sLTE.nv(i));
			std::cout << "level " << i << " from solver: " << s0.nv(i)
			          << ", from LTE:" << sLTE.nv(i) << std::endl;
		}

		// Check if the total density is correct
		double eps = 1e-15;
		CHECK(TemplatedUtils::equalWithinTolerance(s0.nv.sum(), n, eps));
		CHECK(TemplatedUtils::equalWithinTolerance(sLTE.nv.sum(), n, eps));
	}

	SUBCASE("no radiation, low density, should go to ground state")
	{
		double n = 1e-12;
		speciesNv(SpeciesIndex::ine()) = 0;
		speciesNv(SpeciesIndex::inp()) = 0;
		speciesNv(SpeciesIndex::inH()) = 0;
		speciesNv(SpeciesIndex::inH2()) = n;
		const GasStruct gas(T, speciesNv);
		Array specificIntensityv(frequencyv.size());
		Spectrum specificIntensity(frequencyv, specificIntensityv);
		NLevel::Solution s0 = h2l.solveBalance(n, specificIntensity, zerov, zerov, gas);
		EVector nv = s0.nv;
		// Check if some individual levels are indeed 0
		for (size_t i = 1; i < std::min<int>(5, nv.size()); i++)
			WARN_MESSAGE(nv(i) == 0, "nv(" << i << ") = " << nv(i));

		// Check if almost all density is in lowest level
		double eps = 1e-12;
		CHECK(TemplatedUtils::equalWithinTolerance(nv(0), n, eps));
		CHECK((nv.tail(nv.size() - 1).array() < eps).all());
	}
}
