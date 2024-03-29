#include "Array.hpp"
#include "DoctestUtils.hpp"
#include "FreeFree.hpp"
#include "Testing.hpp"
#include "doctest.h"

using namespace RADAGAST;

TEST_CASE("Free-free cooling data vs manual integration")
{
    FreeFree ff;
    Array integrationgridv = Testing::generateGeometricGridv(1000, 1e7, 1e16);
    std::vector<double> Tv = {50., 500., 5000., 50000.};
    for (double T : Tv)
    {
        double ne = 1.e4;
        double np_ne = ne * ne;
        double fromData = ff.cooling(np_ne, T);

        Array gamma_nuv(integrationgridv.size());
        ff.addEmissionCoefficientv(T, integrationgridv, gamma_nuv);

        // emissivity = ne np / 4pi * gamma
        // total emission = 4pi integral(emissivity) = ne np integral(gamma)
        double manuallyIntegrated = np_ne * TemplatedUtils::integrate<double>(integrationgridv, gamma_nuv);

        INFO("cooling data vs manual integration for T = ", T);
        DoctestUtils::checkTolerance("freefree: manually integrated vs from data", manuallyIntegrated, fromData, 1.e-3);
    }
}
