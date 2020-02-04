#include "ChargeDistribution.hpp"
#include "DoctestUtils.hpp"

TEST_CASE("Trivial charge distribution")
{
    ChargeDistribution cd;
    CHECK(cd.zmin() == 0);
    CHECK(cd.zmax() == 0);
    CHECK(cd.value(0) == 1.);
}

TEST_CASE("Charge distribution detailed balance")
{
    ChargeDistribution cd;

    // Some fictional up and down rates. The resulting distribution should satisfy the detailed
    // balance equation: n(i) = n(i - 1) * u(i - 1) / d(i), where n(i), u(i) and d(i) are the
    // density, up rate and down rate of charge i, respectively.
    SUBCASE("flat")
    {
        // When up and down rates are equal and constant, we expect a flat distribution
        auto constant = [](int z) { return 2.; };
        int zmin = -5;
        int zmax = 5;
        cd.calculateDetailedBalance(constant, constant, zmin, zmax);
        int numZ = cd.zmax() - cd.zmin() + 1;
        CHECK(numZ == zmax - zmin + 1);
        double expectedDensity = 1. / numZ;
        for (int z = cd.zmin(); z <= cd.zmax(); z++)
        {
            // This seems to work exactly! Replace by DocTestUtils::checkTolerance if this stop working.
            CHECK(cd.value(z) == expectedDensity);
        }
    }

    SUBCASE("peaked")
    {
        int desiredPeak = 1;
        int zmin = desiredPeak - 3;
        int zmax = desiredPeak + 3;

        // For z < peak, we want the up rate to be higher than the down rate. For z > peak, vice
        // versa. At the peak, the two rates should be roughly (or exactly?) equal.

        double height = 1.;
        double minimum = 0.1;
        double slope = (height - minimum) / (zmax - desiredPeak);

        auto up = [&](int z) { return height - slope * (z - desiredPeak); };
        auto down = [&](int z) { return height + slope * (z - desiredPeak); };
        cd.calculateDetailedBalance(up, down, zmin, zmax);
        for (int z = cd.zmin(); z < desiredPeak; z++) CHECK(cd.value(z) < cd.value(desiredPeak));
        for (int z = desiredPeak + 1; z <= cd.zmax(); z++) CHECK(cd.value(z) < cd.value(desiredPeak));
    }
}
