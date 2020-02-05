#include "DoctestUtils.hpp"
#include "GasInterface.hpp"
#include "GasSolution.hpp"
#include "Testing.hpp"

TEST_CASE("Blackbodies test")
{
    // Test if options don't make it crash
    Array frequencyv = Testing::defaultFrequencyv(100);
    GasModule::GasInterface gi(frequencyv, frequencyv, frequencyv);
    GasModule::GrainInterface gri{};

    std::vector<double> Tv;
    bool warnOnly = false;

    SUBCASE("Very low color temperature - problems are expected here")
    {
        Tv = {66., 112., 225., 450., 900., 1875.};
        warnOnly = true;
    }
    SUBCASE("Normal color temperature - should drive equilibrium")
    {
        Tv = { 3750., 7500., 15000., 30000.};
        warnOnly = false;
    }

    for (double Tc : Tv)
    {
        GasSolution s = gi.solveInitialGuess(1000, Tc, gri);
        double T = s.t();
        CAPTURE(Tc);
        CAPTURE(T);
        DoctestUtils::checkTolerance("equilibrium temp", T, Tc, 0.5, warnOnly);
        WARN_MESSAGE(
            s.nH2() == 0.,
            "without collisional dissociation, H2 won't go to zero (both formation and dissociation likely 0)");
    }
}

TEST_CASE("hff and J V options")
{
    Array frequencyv = Testing::defaultFrequencyv(100);
    GasModule::GasInterface gi = GasModule::GasInterface(frequencyv, frequencyv, frequencyv, "hff2", "8 2");
}

TEST_CASE("hhc option")
{
    Array frequencyv = Testing::defaultFrequencyv(100);
    GasModule::GasInterface gi = GasModule::GasInterface(frequencyv, frequencyv, frequencyv, "hhc");
}
