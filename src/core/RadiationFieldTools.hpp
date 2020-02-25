#ifndef CORE_RADIATIONFIELDTOOLS_HPP
#define CORE_RADIATIONFIELDTOOLS_HPP

#include "Array.hpp"
#include "Constants.hpp"
#include "Spectrum.hpp"

namespace RadiationFieldTools
{
    Array generateBlackbodyv(const Array& frequencyv, double Tc);
    Array generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0);

    // Frequency to wavelength (and the other way around). nu = c / lambda, lambda = c / nu
    Array freqToWavGrid(const Array& frequencyv);
    Array freqToWavSpecificIntensity(const Array& frequencyv, const Array& specificIntensity_nu);

    /** Convert an emissivity (assumed to be in erg / s / cm-3 / Hz) to W / cm-3 / Hz */
    template<typename T> T emissivity_to_SI(T value);

    /** Get the radiation field between 6 and 13.6 eV, in Habing units (1 Habing is 1.6e-3 erg cm-2
        s-1, between 6 and 13.6 eV) */
    double gHabing(const Spectrum& specificIntensity);
    constexpr double nuMinHabing = Constant::LIGHT / (2400 * Constant::ANGSTROM);
    constexpr double nuMaxHabing = Constant::LIGHT / (912 * Constant::ANGSTROM);
}  // namespace RadiationFieldTools

template<typename T> T RadiationFieldTools::emissivity_to_SI(T value)
{
    // erg cm-3 s-1 Hz-1
    //   J  m-3 s-1 Hz-1
    //  -7    6          -> -1
    return 0.1 * value;
}

#endif  // CORE_RADIATIONFIELDTOOLS_HPP
