#ifndef CORE_LEVELSOLVER_HPP
#define CORE_LEVELSOLVER_HPP

#include "EigenAliases.hpp"

namespace LevelSolver
{
    /** Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji +
        B_ji * P_ji + C_ji.

        @note {External processes can influence the level balance by their contribution to the sink
        and source terms. Examples are: level-specific formation (think in the chemical network),
        ionization from specific levels, dissociation from specific levels... so some knowledge
        about the nature of the levels will be necessary. This knowledge is available in for example
        the HModel and BigH2Model classes.

        Sets up F and b using M_ij and the external source terms. Returns the solution as a vector
        [cm-3].

        This is a basic implementation which calls a linear solver from the Eigen library. Different
        types of systems might benefit from other algorithms, based on beforehand knowledge about
        the coefficients. For the H2 model for example, the calculation can be done faster and more
        precisely by using an iterative approach based on the fact that there is no transition data
        between and within the electronically excited levels. */
    EVector statisticalEquilibrium(double totalDensity, const EMatrix& totalTransitionRatesvv, const EVector& sourcev,
                                   const EVector& sinkv, int replaceByConvervationEq = 0);

    /** Solves the statistical equilibrium by looping over the levels one by one, similar to the
        cloudy implementation of H2. This method can be faster than a full matrix inversion when
        there is a large set of levels between which no transitions occur. For H2 for example, the
        electronically excited levels only have transition data with respect to the electronic
        ground state (X) levels, while the X levels also have transition data between them.

        Therefore, this method takes two ranges of indices as extra arguments:

        - fullyConnectedIndices: These indices are connected (with transition coefficients) to each
          other, and to each other level (even the unConnectedIndices).

        - unConnectedIndices: These indices are not connected to each other, and their populations
          only depends on the fullyConnectedIndices.

          In practice, when updating one of the fully connected levels, the populations of all
          levels are included in the calculation. Conversely, when one of the unconnected levels is
          updated, only the fully connected levels' populations is used, as the levels within this
          set are assumed to have no transitions between them. */
    EVector statisticalEquilibrium_iterative(double totalDensity, const EMatrix& totalTransitionRatesvv,
                                             const EVector& sourcev, const EVector& sinkv, const EVector& initialGuessv,
                                             int fullyConnectedCutoff = -1);

    /** Calculate the LTE statistical equilibrium by simply applying the boltzman equations. The
        only thing needed is the energies and the degeneracies of the levels, and a temperature. */
    EVector statisticalEquilibrium_boltzman(double totalDensity, double T, const EVector& ev, const EVector& gv);

} /* namespace LevelSolver */

#endif  // CORE_LEVELSOLVER_HPP
