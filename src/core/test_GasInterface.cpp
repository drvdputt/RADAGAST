#include "DoctestUtils.hpp"
#include "GasInterface.hpp"
#include "GasSolution.hpp"
#include "RadiationFieldTools.hpp"
#include "Testing.hpp"

TEST_CASE("Blackbodies test")
{
    // Test if options don't make it crash
    Array frequencyv = Testing::defaultFrequencyv(300);
    GasModule::GasInterface gi(frequencyv, frequencyv, frequencyv);
    GasModule::GrainInterface gri;
    double nHtotal = 1000;
    std::vector<double> Tv;
    bool warnOnly = false;
    bool withDust = false;

    SUBCASE("Very low color temperature - problems are expected here")
    {
        Tv = {66., 112., 225., 450., 900., 1875.};
        warnOnly = true;
    }
    SUBCASE("Normal color temperature - should drive equilibrium")
    {
        Tv = {3750., 7500., 15000., 30000.};
        warnOnly = false;

        // Try this with and without dust
        SUBCASE("no dust") { withDust = false; }
        SUBCASE("MRNdust")
        {
            withDust = true;
            warnOnly = true;
        }
    }

    for (double Tc : Tv)
    {
        Spectrum specificIntensity(frequencyv, RadiationFieldTools::generateBlackbodyv(frequencyv, Tc));
        if (withDust) Testing::genMRNDust(gri, nHtotal, specificIntensity, true);
        GasSolution s = gi.solveTemperature(nHtotal, specificIntensity, gri);
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
