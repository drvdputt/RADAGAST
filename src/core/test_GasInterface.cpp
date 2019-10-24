#include "DoctestUtils.hpp"
#include "GasInterface.hpp"
#include "GasSolution.hpp"
#include "Testing.hpp"

TEST_CASE("Gas interface")
{
	// Test if options don't make it crash
	Array frequencyv = Testing::defaultFrequencyv(100);
	std::unique_ptr<GasModule::GasInterface> gip{nullptr};

	SUBCASE("Full model blackbody test")
	{
		std::vector<double> Tv = {1875., 3750., 7500., 15000., 30000.};
		for (double Tc : Tv)
		{
			GasModule::GasInterface gi(frequencyv, frequencyv, frequencyv);
			GasModule::GrainInterface gri{};
			GasSolution s = gi.solveInitialGuess(1000, Tc, gri);
			double T = s.t();
			CAPTURE(Tc);
			CAPTURE(T);
			DoctestUtils::checkTolerance("equilibrium temp vs color temp", T, Tc, 0.5);
			WARN_MESSAGE(s.nH2() == 0, "without collisional dissociation, H2 won't go to "
				     "zero (both formation and dissociation likely 0)");
		}
	}
	SUBCASE("hhc option")
	{
		GasModule::GasInterface gi = GasModule::GasInterface(frequencyv, frequencyv,
		                                                     frequencyv, "hhc");
	}
	SUBCASE("hff and J V options")
	{
		GasModule::GasInterface gi = GasModule::GasInterface(frequencyv, frequencyv,
		                                                     frequencyv, "hff2", "8 2");
	}
}
