#include "RadiationFieldTools.hpp"
#include "SpecialFunctions.hpp"
#include <vector>

Array RadiationFieldTools::generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0)
{
    Array I_nu(frequencyv.size());
    for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
        I_nu[iFreq] = SpecialFunctions::planck(frequencyv[iFreq], Tc);

    // Rescale to the desired G0
    double currentG0 = gHabing(Spectrum(frequencyv, I_nu));
    I_nu *= G0 / currentG0;
    return I_nu;
}

Array RadiationFieldTools::freqToWavGrid(const Array& frequencyv)
{
    size_t numWav = frequencyv.size();
    Array wavelengthv(numWav);
    for (size_t iWav = 0; iWav < numWav; iWav++) wavelengthv[iWav] = Constant::LIGHT / frequencyv[numWav - 1 - iWav];
    return wavelengthv;
}

Array RadiationFieldTools::freqToWavSpecificIntensity(const Array& frequencyv, const Array& specificIntensity_nu)
{
    Array I_lambda(frequencyv.size());
    for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
        I_lambda[I_lambda.size() - iFreq - 1] =
            specificIntensity_nu[iFreq] * frequencyv[iFreq] * frequencyv[iFreq] / Constant::LIGHT;
    return I_lambda;
}

double RadiationFieldTools::gHabing(const Spectrum& specificIntensity)
{
    double intensity = specificIntensity.average(nuMinHabing, nuMaxHabing) * (nuMaxHabing - nuMinHabing);
    return Constant::FPI * intensity / Constant::HABING_FLUX;
}
