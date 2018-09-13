#include "LevelSolver.h"

EVector LevelSolver::statisticalEquilibrium(double totalDensity,
                                            const EMatrix& totalTransitionRatesvv,
                                            const EVector& sourcev, const EVector& sinkv)
{
	// For an explanation of the algorithm, see gasphysics document

	// The transition matrix T_ij gives on on each row i, the departure rate from the level
	// i towards each other level. Each column j thus represents the arrival rates into
	// level j, from each other level. Multiplying with a population row vector n_i on the
	// left, will give a new vector containing the the total arrival rates in each level. We
	// work with the transpose here, so that Mvv * nv gives a column vector containing the
	// total arrival rates.
	EMatrix Mvv = totalTransitionRatesvv.transpose();

	// Take the sum of the transition coefficients from each level i, resulting in the
	// total fraction of this level that transitions into other levels per unit time.
	EMatrix departureDiagonal = Mvv.colwise().sum().asDiagonal();
	Mvv -= departureDiagonal; // Fij
	Mvv -= sinkv.asDiagonal(); // Fij - diag[d_i]
	EVector f(-sourcev); // -f_i

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = EVector::Ones(Mvv.cols());
	f(chooseConsvEq) = totalDensity;

#ifdef PRINT_LEVEL_MATRICES
	DEBUG("System to solve:\n" << Mvv << " * nv\n=\n" << f << endl << endl);
#endif
	// Call the linear solver for the system sum_j (Fij - diag[d_i]) * n_j = -f_i
	EVector nv = Mvv.colPivHouseholderQr().solve(f);

	// Put populations = 0 if they were negative due to precision issues
	nv = nv.array().max(0);
	return nv;
}

EVector LevelSolver::statisticalEquilibrium_iterative(double totalDensity,
                                                      const EMatrix& totalTransitionRatesvv,
                                                      const EVector& sourcev,
                                                      const EVector& sinkv,
                                                      const EVector& initialGuessv)
{
	EMatrix Mvv = totalTransitionRatesvv.transpose();
	EVector fracDestructionRatev = Mvv.colwise().sum().transpose() + sinkv;
}
