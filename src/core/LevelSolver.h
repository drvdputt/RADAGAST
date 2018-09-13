#ifndef GASMODULE_GIT_SRC_LEVELSOLVER_H_
#define GASMODULE_GIT_SRC_LEVELSOLVER_H_

#include NLevel

namespace LevelSolver
{
/** Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji +
    B_ji * P_ji + C_ji. 
    
    @note {External processes can influence the level balance by their contribution to the sink
    and source terms. Examples are: level-specific formation (think in the chemical network),
    ionization from specific levels, dissociation from specific levels... so some knowledge
    about the nature of the levels will be necessary. Subclasses of NLevel can provide an
    interface to whatever that knowledge may be (usually a mapping from quantum numbers to level
    index), but the general NLevel implementation is completely oblivious to this, and hence
    provides no source and sink coefficients.}

    Sets up F and b using M_ij and the external source terms.
    Returns the solution as a vector [cm-3].

    This is a basic implementation which calls a linear solver from the Eigen library. Different
    types of systems might benefit from other algorithms, based on beforehand knowledge about
    the coefficients. For the H2 model for example, the calculation can be done faster and more
    precisely by using an iterative approach based on the fact that there is no transition data
    between and within the electronically excited levels. */
EVector statisticalEquilibrium(double totalDensity, const EMatrix& totalTransitionRatesvv,
                               const EVector& sourcev, const EVector& sinkv);

EVector statisticalEquilibrium_iterative(double totalDensity,
                                         const EMatrix& totalTransitionRatesvv,
                                         const EVector& sourcev, const EVector& sinkv, const EVector& initialGuessv);

} /* namespace LevelSolver */

#endif /* GASMODULE_GIT_SRC_LEVELSOLVER_H_ */
