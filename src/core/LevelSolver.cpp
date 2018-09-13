#include "LevelSolver.h"
#include "DebugMacros.h"
#include "Error.h"

EVector LevelSolver::statisticalEquilibrium(double totalDensity,
                                            const EMatrix& totalTransitionRatesvv,
                                            const EVector& sourcev, const EVector& sinkv,
                                            int replaceByConservationEq)
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
	Mvv.row(replaceByConservationEq) = EVector::Ones(Mvv.cols());
	f(replaceByConservationEq) = totalDensity;

#ifdef PRINT_LEVEL_MATRICES
	DEBUG("System to solve:\n" << Mvv << " * nv\n=\n" << f << std::endl << std::endl);
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
                                                      const EVector& initialGuessv,
                                                      int fullyConnectedCutoff)
{
	EVector nv = initialGuessv;
	int numLv = initialGuessv.size();
	EMatrix Mvv = totalTransitionRatesvv.transpose();

	// Fractional destruction rate (in s-1) stays constant when populations are adjusted
	// (sum over the row indices by doing a colwise reduction. Need to transpose because
	// summing each column gives a row vector, while sinkv is column one.
	EVector fracDestructionRatev = Mvv.colwise().sum().transpose() + sinkv;

	// The comments here are written in context of H2, since that is the main reason this
	// algorithm exists

	// Get the indices that cover the electronic ground state. If this information is not
	// given, treat all levels as fully connected.
	int Xstart = 0;
	int Xstop = fullyConnectedCutoff != -1 ? fullyConnectedCutoff : numLv;
	int Estart = Xstop + 1;
	int Estop = numLv;

	// The algorithm apparently works better when the iterations happen from high to low
	// energies, so we go backwards in our loops.

	// Iterate until converged
	double maxDeltaX = 1.e-4 * totalDensity;
	double maxDeltaAll = 1.e-6 * totalDensity;
	double maxFracDelta = 1.e-6;
	size_t counter{0};
	const int max_iterations = 2001;
	while (counter < max_iterations)
	{
		counter++;

		// The previous nv, for convergence checking
		EVector previousNv = nv;
		EVector deltav = EVector::Zero(nv.size());

		// Sweep over ground state
		for (int i = Xstart; i < Xstop; i++)
		{
			// Sum Mij nj, with j running over all other levels. TODO: consider
			// making the assumption that the X states are contigous in the eigen
			// index space. In that case, I can simply do a matrix operation here,
			// over a slice, leading to possible optimization.
			double creationRate = sourcev(i) + Mvv.row(i) * nv;
			nv(i) = creationRate <= 0 ? 0 : creationRate / fracDestructionRatev(i);
		}
		/* Renormalize because the algorithm has no sum rule, */
		nv *= totalDensity / nv.sum();

		// We will cut this iteration short as long as the ground state has not
		// converged.
		bool loopBack = false;
		for (int i = Xstart; i < Xstop; i++)
		{
			// TODO: vectorize this (maybe work with a contiguous 'range' of X
			// indices, instead of a list
			deltav(i) = abs(nv(i) - previousNv(i));
			if (!loopBack && deltav(i) > maxDeltaX)
			{
				// DEBUG("index " << i << " delta " << deltav(i) / maxDeltaX << std::endl);
				// Do not break here. We might need the whole deltav further down
				loopBack = true;
			}
		}
		if (loopBack)
		{
			if (!(counter % 1000))
			{
				DEBUG("Solving h2... " << counter << "iterations" << std::endl);
				// DEBUG("h2Levelv = \n" << nv << std::endl);
			}
			continue;
		}

		// If the ground state has more or less converged, we will also start sweeping
		// over the other states
		for (int i = Estart; i < Estop; i++)
		{
			double creationRate = sourcev(i);
			// Only sum over the ground state here
			for (int j = Xstart; j < Xstop; j++)
				creationRate += Mvv(i, j) * nv(j);
			
			nv(i) = creationRate <= 0 ? 0 : creationRate / fracDestructionRatev(i);
		}
		nv *= totalDensity / nv.sum();

		// Overall convergence check

		// Absolute change
		for (int i = Estart; i < Estop; i++)
			// Again, this would best be vectorized
			deltav(i) = abs(nv(i) - previousNv(i));

		bool thresconv = (deltav.array() < maxDeltaAll).all();

		// Relative change
		bool fracconv = true;
		for (int i = 0; i < nv.size() && fracconv; i++)
		{
			double df = 0;
			// finite/0 means not converged
			if (previousNv(i) <= 0 && nv(i) > 0)
				df = 2 * maxFracDelta;
			else
				df = abs(nv(i) / previousNv(i) - 1.);

			if (df > maxFracDelta)
				fracconv = false;
		}

		if (!(counter % 1000))
		{
			DEBUG("Solving h2... " << counter << "iterations" << std::endl);
			DEBUG("thresconv = " << thresconv << " fracconv = " << fracconv
			                     << std::endl);
		}

		bool converged = thresconv && fracconv;
		if (converged)
			break;
	}
	if (nv.array().isNaN().any())
		Error::runtime("nan in H2 level solution");
	DEBUG("Solved H2 in " << counter << " iterations" << std::endl);
	// DEBUG("h2Levelv = \n" << nv << std::endl);
	// EVector explicitNv = NLevel::solveRateEquations(n, BPvv, Cvv, sourcev, sinkv,
	// chooseConsvEq, gas);
	return nv;
}
