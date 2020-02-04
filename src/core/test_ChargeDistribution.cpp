#include "ChargeDistribution.hpp"
#include "DoctestUtils.hpp"

TEST_CASE("Trivial charge distribution")
{
    ChargeDistribution cd;
    CHECK(cd.zmin() == 0);
    CHECK(cd.zmax() == 0);
    CHECK(cd.value(0) == 1.);
}

TEST_CASE("Flat charge distribution")
{
    ChargeDistribution cd;

    // Some fictional up and down rates. The resulting distribution should satisfy the detailed
    // balance equation: n(i) = n(i - 1) * u(i - 1) / d(i), where n(i), u(i) and d(i) are the
    // density, up rate and down rate of charge i, respectively.
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

TEST_CASE("Peaked charge distribution")
{
    ChargeDistribution cd;
    int desiredPeak, zmin, zmax, maxCharges;
    double height, minimum, slope;

    SUBCASE("simple")
    {
        desiredPeak = 1;
        zmin = desiredPeak - 3;
        zmax = desiredPeak + 3;
        maxCharges = 0;
        // For z < peak, we want the up rate to be higher than the down rate. For z > peak, vice
        // versa. At the peak, the two rates should be roughly (or exactly?) equal.
        height = 1.;
        minimum = 0.1;
    }

    SUBCASE("wide and strong - trimming needed")
    {
        desiredPeak = -300;
        zmin = -600;
        zmax = 100;
        maxCharges = 0;
        height = 100.;
        minimum = 0.1;
    }

    SUBCASE("wide and strong - enforce numCharges")
    {
        desiredPeak = -300;
        zmin = -600;
        zmax = 100;
        maxCharges = 50;
        height = 100.;
        minimum = 0.1;
    }

    // Linear functions equal to minimum at zmin (for down rate) or zmax (for up rate). At
    // desiredPeak, they cross.
    slope = (height - minimum) / (zmax - desiredPeak);
    auto up = [&](int z) { return height - slope * (z - desiredPeak); };
    auto down = [&](int z) { return height + slope * (z - desiredPeak); };

    cd.calculateDetailedBalance(up, down, zmin, zmax);

    // Check that the peak is where expected
    for (int z = cd.zmin(); z < desiredPeak; z++) CHECK(cd.value(z) < cd.value(desiredPeak));
    for (int z = desiredPeak + 1; z <= cd.zmax(); z++) CHECK(cd.value(z) < cd.value(desiredPeak));

    // If automatic trimming has happened, check that ratio of the density at the edge of the
    // distribution is small enough compared to the densityt at the peak. Currently, the ratio
    // where trimming occurs is hardcoded. Check the source of calculateDetailedBalance. Maybe make
    // this fraction an option later.
    if (cd.zmin() != zmin)
    {
        CHECK(cd.zmin() > zmin);
        CHECK(cd.value(cd.zmin()) / cd.value(desiredPeak) < 0.1);
    }
    if (cd.zmax() != zmax)
    {
        CHECK(cd.zmax() < zmax);
        CHECK(cd.value(cd.zmax()) / cd.value(desiredPeak) < 0.1);
    }
}
