#include "doctest.h"

#include "GasStruct.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "SpeciesIndex.h"
#include "Testing.h"

TEST_CASE("H2-specific algorithm")
{
	H2FromFiles hff(2, 2);
	H2Levels h2l(std::make_shared<const H2FromFiles>(hff));

	double T = 100;
	double n = 100;
	EVector speciesNv = EVector::Zero(SpeciesIndex::size());
	speciesNv(SpeciesIndex::ine()) = 1;
	speciesNv(SpeciesIndex::inp()) = 1;
	speciesNv(SpeciesIndex::inH()) = 1;
	speciesNv(SpeciesIndex::inH2()) = n;
	const GasStruct gas(T, speciesNv);

	EVector zerov = EVector::Zero(h2l.numLv());

	Array frequencyv = Testing::generateGeometricGridv();

	double epsFrac = 1e-12;

	SUBCASE("no radiation")
	{
		Array specificIntensityv(frequencyv.size());
		Spectrum specificIntensity(frequencyv, specificIntensityv);
		NLevel::Solution s0 = h2l.solveBalance(n, specificIntensity, zerov, zerov, gas);
		NLevel::Solution sLTE = h2l.solveLTE(n, specificIntensity, gas);

		for (size_t i = 0; i < s0.nv.size(); i++)
			CHECK_MESSAGE(TemplatedUtils::equalWithinTolerance(s0.nv(i), sLTE.nv(i),
			                                                   epsFrac),
			              "level " << i << " from solver: " << s0.nv(i)
				      << ", from LTE:" << sLTE.nv(i));
	}
}
