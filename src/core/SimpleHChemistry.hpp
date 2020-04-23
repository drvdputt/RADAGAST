#ifndef CORE_SIMPLEHCHEMISTRY_HPP
#define CORE_SIMPLEHCHEMISTRY_HPP

#include "Chemistry.hpp"

namespace GasModule
{
    class Spectrum;

    class SimpleHChemistry : public Chemistry
    {
    public:
        /** Sets up a chemical network describing only the most basic reactions between e-, H+, H
            and H2.

            - photoionization of H [s-1]

            - ionization by collission with electron [cm3 s-1]

            - radiative recombination [cm3 s-1]

            - dissociation of H2 after excitation [s-1]

            - H2 formation (twice the H2 formation rate, as we use 1 H -> 0.5 H2 to make the
              reaction scale linearly with nH) [s-1] */
        SimpleHChemistry();

        /** Calculate the rate coefficients for each reaction. Multiplying with the right density
            factors of the reactants will give a total rate [cm-3 s-1]. The unit of these
            coefficients is hence [cm(-3+n) s-1] where n is usually the number of reactant
            particles. Since the calculation of the H2 formation and destruction rate are
            complicated, the are expected to be calculated somewhere else, and passed as arguments
            here. This function simply fills them in into the right spot of the k-vector. */
        EVector rateCoeffv(double T, const Spectrum& meanIntensity, double kDissFromH2Levels,
                           double kH2FormationGrain) const;
    };
}
#endif  // CORE_SIMPLEHCHEMISTRY_HPP
