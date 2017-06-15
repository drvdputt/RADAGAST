#include "UnitTests.h"
#include "Constants.h"
#include "HydrogenLevels.h"
#include "Testing.h"

#include <vector>

using namespace std;

void UTest::testHardcodedHydrogen()
{
	auto freqv = Testing::generateGeometricGridv(150, Constant::LIGHT / (1e4 * Constant::UM_CM),
	                                             Constant::LIGHT / (0.01 * Constant::UM_CM));
	double Tc = 20000;
	double G0 = 1e0;
	double n = 1e1;
	HydrogenLevels hl(Array(freqv.data(), freqv.size()));
	Array specificIntensityv = Testing::generateSpecificIntensityv(
	                vector<double>(begin(freqv), end(freqv)), Tc, G0);

	// TODO: write an actual test
}

