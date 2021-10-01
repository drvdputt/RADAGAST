#include "DoctestUtils.hpp"
#include "Ionization.hpp"
#include "RadiationFieldTools.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"
#include "Ionization.hpp"

using namespace RADAGAST;

TEST_CASE("Ionization threshold not in grid")
{
    // There was a bug where the code segfaulted when a wavelength grid was chosen for which the
    // H ionization threshold frequency was not in its bounds.

    Array frequencyv;
    double nu0 = Ionization::THRESHOLD;
    bool shouldBeZero;

    // two subcases we want to test:
    SUBCASE("max frequency < threshold")
    {
        frequencyv = Testing::generateGeometricGridv(200, nu0 / 4, nu0 / 1.1);
        shouldBeZero = true;
    }
    SUBCASE("threshold < min frequency")
    {
        frequencyv = Testing::generateGeometricGridv(200, nu0 * 1.1, nu0 * 2);
        shouldBeZero = false;
    }

    const Array& meanIntensityv = RadiationFieldTools::generateSpecificIntensityv(frequencyv, 10000, 100);
    Spectrum meanIntensity(frequencyv, meanIntensityv);
    double photoRate = Ionization::photoRateCoeff(meanIntensity);

    if (shouldBeZero)
        CHECK(photoRate == 0.);
    else
        CHECK(photoRate > 0.);
}
