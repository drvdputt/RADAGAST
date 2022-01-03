#include "DoctestUtils.hpp"
#include "GasInterface.hpp"
#include "GasSolution.hpp"
#include "GasState.hpp"
#include "RadiationFieldTools.hpp"
#include "Testing.hpp"

using namespace RADAGAST;

TEST_CASE("Blackbodies test")
{
    // Test if options don't make it crash
    Array frequencyv = Testing::defaultFrequencyv(300);
    RADAGAST::GasInterface gi(frequencyv, frequencyv, frequencyv);
    double nHtotal = 1000;
    std::vector<double> Tv;
    bool warnOnly = false;
    bool withDust = false;

    SUBCASE("Very low color temperature - problems are expected here")
    {
        Tv = {100., 300., 900., 1875.};
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
            // Note: technically, the maximum grain temperature should be increased for this
            // test (ideally, it ends up being == Tc).
            withDust = true;
            warnOnly = true;
        }
    }

    for (double Tc : Tv)
    {
        Spectrum meanIntensity(frequencyv, RadiationFieldTools::generateBlackbodyv(frequencyv, Tc));
        RADAGAST::GrainInterface gri;
        if (withDust) Testing::genMRNDust(gri, nHtotal, meanIntensity, true);
        GasSolution s = gi.solveTemperature(nHtotal, meanIntensity, 1., gri);
        double T = s.t();
        CAPTURE(Tc);
        CAPTURE(T);
        // Ideally, the equilibrium temperature should approach the color temperature. Note that
        // this will only work if all the implemented processes have their thermodynamic
        // counterpart implemented too, which is probably not the case for the dust
        // photoelectric/collision processes. Therefore, warnOnly is turned on when we have dust
        // present. It works surprisingly well when running without dust though.
        DoctestUtils::checkTolerance("equilibrium temp", T, Tc, 0.5, warnOnly);
        if (!withDust)
        {
            // Without dust, we expect no H2 to be present
            DoctestUtils::checkRange("zero h2 dens", s.nH2(), 0., 1.e-17 * nHtotal, true);
        }
    }
}

TEST_CASE("zero radiation field")
{
    Array frequencyv = Testing::defaultFrequencyv(300);
    Array meanIntensityv(frequencyv.size());
    Spectrum meanIntensity(frequencyv, meanIntensityv);
    RADAGAST::GasInterface gi(frequencyv, frequencyv, frequencyv);
    double nHtotal = 0.;
    RADAGAST::GrainInterface gri;

    SUBCASE("zero dens") {}
    SUBCASE("nonzero dens")
    {
        nHtotal = 100.;
        SUBCASE("zero grains") {}
        SUBCASE("nonzero grains") { Testing::genMRNDust(gri, nHtotal, meanIntensity, true); }
    }

    GasSolution s = gi.solveTemperature(nHtotal, meanIntensity, 1., gri);
    // no checks for now, just make sure it doesn't crash
}

TEST_CASE("serialization")
{
    // need to make sure this test fails if the current assumptions for the very simple
    // serialize implementation are violated.
    auto gi = Testing::genFullModel();
    GrainInterface gri;
    GasState gs;
    auto iv = RadiationFieldTools::generateSpecificIntensityv(gi.iFrequencyv(), 1e4, 100);
    gi.updateGasState(gs, 1000, iv, 1., gri);

    auto blob = gi.serialize(gs);
    // hardcode this, to warn me later, when I change the chemistry (number of species), or the
    // implementation of GasStae. With the way it's implemented now, this should be true
    CHECK(blob.size() == 6);
    CHECK(gs.density(0) == blob[1]);

    // manually edit blob to test deserialize, and check again
    blob[1] = 0;
    gi.deserialize(gs, blob);
    CHECK(gs.density(0) == blob[1]);

    // check info function
    auto info = gi.serializationInfo();
    CHECK(info.size() == blob.size());
}
