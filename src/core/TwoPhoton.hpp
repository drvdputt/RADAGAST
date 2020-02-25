#ifndef CORE_TWOPHOTON_HPP
#define CORE_TWOPHOTON_HPP

#include "Array.hpp"

namespace TwoPhoton
{
    /** Calculates the two-photon continuum for H on a grid of frequencies, given the population
        (in cm-3) of the n2s level of H. This function is defined separately, since n2s is the only
        parameter that is needed. [erg s-1 cm-3 Hz-1]. */
    Array emissivityv(const Array& eFrequencyv, double n2s);
}

#endif  // CORE_TWOPHOTON_HPP
