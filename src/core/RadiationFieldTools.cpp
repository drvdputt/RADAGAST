#include "RadiationFieldTools.hpp"
#include "Functions.hpp"
#include <vector>

namespace RADAGAST
{
    Array RadiationFieldTools::generateBlackbodyv(const Array& frequencyv, double Tc)
    {
        Array I_nu(frequencyv.size());
        for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
            I_nu[iFreq] = Functions::planck(frequencyv[iFreq], Tc);
        return I_nu;
    }

    Array RadiationFieldTools::generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0)
    {
        const Array& I_nu = generateBlackbodyv(frequencyv, Tc);
        double currentG0 = gHabing(Spectrum(frequencyv, I_nu));
        return I_nu * G0 / currentG0;
    }

    Array RadiationFieldTools::freqToWavGrid(const Array& frequencyv)
    {
        size_t numWav = frequencyv.size();
        Array wavelengthv(numWav);
        for (size_t iWav = 0; iWav < numWav; iWav++)
            wavelengthv[iWav] = Constant::LIGHT / frequencyv[numWav - 1 - iWav];
        return wavelengthv;
    }

    Array RadiationFieldTools::freqToWavSpecificIntensity(const Array& frequencyv, const Array& meanIntensity_nu)
    {
        Array I_lambda(frequencyv.size());
        for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
            I_lambda[I_lambda.size() - iFreq - 1] =
                meanIntensity_nu[iFreq] * frequencyv[iFreq] * frequencyv[iFreq] / Constant::LIGHT;
        return I_lambda;
    }

    double RadiationFieldTools::gHabing(const Spectrum& meanIntensity)
    {
        double intensity = meanIntensity.average(nuMinHabing, nuMaxHabing) * (nuMaxHabing - nuMinHabing);
        return Constant::FPI * intensity / Constant::HABING_FLUX;
    }
}
