#include "DoctestUtils.h"
#include "GasStruct.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "LevelSolver.h"
#include "SpeciesIndex.h"
#include "Testing.h"

TEST_CASE("H2-specific algorithm")
{
	H2FromFiles hff(15, 3);
	EVector ev = hff.ev();
	H2Levels h2l(std::make_shared<const H2FromFiles>(hff));

	double T = 500;

	EVector speciesNv = EVector::Zero(SpeciesIndex::size());
	EVector zerov = EVector::Zero(h2l.numLv());
	Array frequencyv = Testing::generateGeometricGridv();

	double epsFrac = 1e-2;

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

		NLevel::Solution s0 = h2l.customSolution(n, gas, specificIntensity);
		NLevel::Solution sLTE = h2l.solveLTE(n, gas);

		// This test seems to work reasonably for the first three levels
		for (size_t i = 0; i < std::min<int>(3, s0.nv.size()); i++)
			DoctestUtils::checkTolerance("level " + std::to_string(i) +
			                                             " from solver vs LTE",
			                             s0.nv(i), sLTE.nv(i), epsFrac);

		// Check if the total density is correct
		double eps = 1e-15;
		DoctestUtils::checkTolerance("sum of s0.nv", s0.nv.sum(), n, eps);
		DoctestUtils::checkTolerance("sum of sLTE.nv", sLTE.nv.sum(), n, eps);
	}

	SUBCASE("no radiation, low density, should go to ground state")
	{
		double n = 1e-12;
		double eps = 1e-12;

		speciesNv(SpeciesIndex::ine()) = 0;
		speciesNv(SpeciesIndex::inp()) = 0;
		speciesNv(SpeciesIndex::inH()) = 0;
		speciesNv(SpeciesIndex::inH2()) = n;
		const GasStruct gas(T, speciesNv);
		Array specificIntensityv(frequencyv.size());
		Spectrum specificIntensity(frequencyv, specificIntensityv);
		NLevel::Solution s0 = h2l.customSolution(n, gas, specificIntensity);
		EVector nv = s0.nv;

		// Check if some individual levels are indeed close to 0
		for (size_t i = 2; i < std::min<int>(5, nv.size()); i++)
			DoctestUtils::checkRange("level " + std::to_string(i), nv(i), 0,
			                         n * eps);

		// Check if almost all density is in lowest (two) level(s) (ortho+para)
		DoctestUtils::checkTolerance("ground level(s)", nv(0) + nv(1), n, eps);

		CAPTURE("ground level (energy "
		        << ev(0) << ") " << nv(0) << " \n level 1 (energy " << ev(1) << ") "
		        << nv(1) << " \n level 3 (energy " << ev(3) << ") " << nv(3) << "\n");

		// Check if all but the lowest two level populations are small compared to total
		CHECK((nv.tail(nv.size() - 2).array() < eps * n).all());
	}

	// TODO: This is not working as expected. Of course, the above examples do not depend
	// on radiation field. I need a way to check that the effect on radiation on level
	// populations is hanled correctly.
	SUBCASE("Blackbody radiation, no collisions, should also go to LTE?")
	{
		double n = 1;
		double T = 250;
		speciesNv(SpeciesIndex::ine()) = 0;
		speciesNv(SpeciesIndex::inp()) = 0;
		speciesNv(SpeciesIndex::inH()) = 0;
		speciesNv(SpeciesIndex::inH2()) = n;
		const GasStruct gas(T, speciesNv);
		Array specificIntensityv(frequencyv.size());
		for (size_t i = 0; i < frequencyv.size(); i++)
			specificIntensityv[i] = SpecialFunctions::planck(frequencyv[i], T);
		Spectrum specificIntensity(frequencyv, specificIntensityv);

		NLevel::Solution s0 = h2l.customSolution(n, gas, specificIntensity);
		NLevel::Solution sLTE = h2l.solveLTE(n, gas);
		for (size_t i = 0; i < std::min<int>(4, s0.nv.size()); i++)
			DoctestUtils::checkTolerance("level pop vs LTE", s0.nv(i), sLTE.nv(i),
			                             epsFrac, true);

		// Check if the total density is correct
		double eps = 1e-15;
		DoctestUtils::checkTolerance("s0.nv.sum() vs n", s0.nv.sum(), n, eps);
		DoctestUtils::checkTolerance("sLTE.nv.sum() vs n", sLTE.nv.sum(), n, eps);
	}
}
