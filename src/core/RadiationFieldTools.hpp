#ifndef CORE_RADIATIONFIELDTOOLS_HPP
#define CORE_RADIATIONFIELDTOOLS_HPP

#include "Array.hpp"
#include "Constants.hpp"
#include "Spectrum.hpp"

namespace GasModule
{
    namespace RadiationFieldTools
    {
        /** Evaluate blackbody intensity on the given grid [erg s-1 Hz-1 cm-2 sr-1]*/
        Array generateBlackbodyv(const Array& frequencyv, double Tc);

        /** Evaluate blackbody, and renormalize it to correspond to a certain G0. */
        Array generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0);

        /** Convert list of frequencies to wavelengths (and the other way around). nu = c / lambda,
            lambda = c / nu */
        Array freqToWavGrid(const Array& frequencyv);

        /** Convert specific intensity from per frequency to per wavelength units */
        Array freqToWavSpecificIntensity(const Array& frequencyv, const Array& meanIntensity_nu);

        /** Convert an emissivity (assumed to be in erg / s / cm-3 / Hz) to W / cm-3 / Hz */
        template<typename T> T emissivity_to_SI(T value);

        /** Get the radiation field between 6 and 13.6 eV (2066 and 912 Angstrom), relative to
            the Habing unit radiation field (see @c Constants) */
        double gHabing(const Spectrum& meanIntensity);
        constexpr double nuMinHabing = 6. * Constant::EV / Constant::PLANCK;
        constexpr double nuMaxHabing = Constant::LIGHT / (912. * Constant::ANGSTROM);
    }  // namespace RadiationFieldTools

    template<typename T> T RadiationFieldTools::emissivity_to_SI(T value)
    {
        // erg cm-3 s-1 Hz-1
        //   J  m-3 s-1 Hz-1
        //  -7    6          -> -1
        return 0.1 * value;
    }
}
#endif  // CORE_RADIATIONFIELDTOOLS_HPP
