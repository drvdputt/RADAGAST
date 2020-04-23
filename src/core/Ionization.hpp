#ifndef CORE_IONIZATION_HPP
#define CORE_IONIZATION_HPP

#include "Spectrum.hpp"

namespace GasModule
{
    namespace Ionization
    {
        /** Solves the ionization balance (in the nebular approximation, i.e. this function assumes
            that all hydrogen is in the ground state). */
        double solveBalance(double nH, double T, const Spectrum& meanIntensity);

        /** Photoionization rate coefficient. Integrates over the spectrum. Multiply with neutral
            density to get total rate. */
        double photoRateCoeff(const Spectrum& meanIntensity);

        /** Photoionization cross section in cm2 */
        double crossSection(double frequency);

        /** Collisional ionization rate coefficient from the ground state in cm3 / s */
        double collisionalRateCoeff(double T);

        /** Total radiative recombination rate coefficient in cpm3 / s */
        double recombinationRateCoeff(double T);

        /** Heating due to thermalization of freed electrons (erg s-1). Multiply with nH to get
            heating rate. */
        double heatingPerH(const Spectrum& meanIntensity);

        /** Kinetic energy lost during recombination (basically the photon energy, minus the
            binding energy contribution. */
        double cooling(double nH, double np, double ne, double T);

        const double THRESHOLD = 3.28984196e15;  // Hertz;
    }
}
#endif  // CORE_IONIZATION_HPP
