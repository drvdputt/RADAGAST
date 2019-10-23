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
		GasModule::GasInterface gi(frequencyv, frequencyv, frequencyv);
		GasModule::GrainInterface gri{};
		GasSolution s = gi.solveInitialGuess(1000, 5000, gri);
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
