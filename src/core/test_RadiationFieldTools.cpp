#include "DoctestUtils.hpp"
#include "RadiationFieldTools.hpp"
#include "Testing.hpp"

TEST_CASE("Calculation of G")
{
	Array frequencyv = Testing::defaultFrequencyv(2000);

	double Tc;
	SUBCASE("Tc = 30000") { Tc = 30000; }
	SUBCASE("Tc = 15000") { Tc = 15000; }
	SUBCASE("Tc = 7500") { Tc = 7500; }

	double G0 = 1;
	Array specificIntensityv =
	                RadiationFieldTools::generateSpecificIntensityv(frequencyv, Tc, G0);

	Spectrum specificIntensity(frequencyv, specificIntensityv);
	double G = RadiationFieldTools::gHabing(specificIntensity);
	DoctestUtils::checkTolerance("G vs G0", G, G0, 1e-2);
}
