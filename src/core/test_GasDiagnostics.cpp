#include "DoctestUtils.hpp"
#include "GasDiagnostics.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "RadiationFieldTools.hpp"
#include "Testing.hpp"

using namespace RADAGAST;

TEST_CASE("GasDiagnostics maybe crashes")
{
    // This is mostly a test to make sure that it doesn't crash under different configurations,
    // since the diagnostics will be digging into specific parts of the simulations. Originally
    // written because of suspected crash with simple H2 model.
    GasInterface gi = Testing::genFullModel();
    GasState gs;
    GasDiagnostics gd;
    GrainInterface gri;
    Array inu = RadiationFieldTools::generateSpecificIntensityv(gi.iFrequencyv(), 10000, 100);

    // Ideally, I'd want this test to run with both the simple and the big H2 model. But
    // currently, that option is set at compile time. But this test should at least help me
    // debug things.
    gi.updateGasState(gs, 10000, inu, 1., gri, &gd);
}
