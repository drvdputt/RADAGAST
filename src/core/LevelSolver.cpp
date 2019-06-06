#include "LevelSolver.h"
#include "Constants.h"
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
	// Tij is the transition rate from i to j. Mij = Tji means the arrival in state i from
	// j. This form is handy because then you can do a left-multiplication with the nv
	// column vector to arrive at the total increase rate for each level. Storing Mvv row
	// major gives some extra speed, because the rows are then contiguous in memory, which
	// helps when multiplying with nv.
	EMatrixRM Mvv = totalTransitionRatesvv.transpose();

	// Fractional destruction rate (in s-1) stays constant when populations are adjusted
	// (sum over the row indices by doing a colwise reduction. Need to transpose because
	// summing each column gives a row vector, while sinkv is column one.
	EVector fracDestructionRatev = totalTransitionRatesvv.rowwise().sum() + sinkv;

	// The comments here are written in context of H2, since that is the main reason this
	// algorithm exists

	// Get the indices that cover the electronic ground state. If this information is not
	// given, treat all levels as fully connected.
	int startX = 0;
	int stopX = fullyConnectedCutoff != -1 ? fullyConnectedCutoff : numLv;
	int numX = stopX - startX;
	int startE = stopX;
	int stopE = numLv;
	int numE = stopE - startE;

	// Convergence criteria. For the ground state levels, we demand that a certain absolute
	// precision is reached. For the excited levels, a fractional convergence criterium
	// suffices.
	double maxDeltaX = 1.e-3 * totalDensity;
	double maxFracDelta = 1.e-3;

	bool converged = false;
	size_t counter{0};
	const int max_iterations = 2001;
	while (!converged && counter < max_iterations)
	{
		counter++;

		EVector previousNv = nv;

		// Sweep over ground state (this can also be done with a block operation, but
		// this is faster for some reason)
		for (int i = startX; i < stopX; i++)
		{
			// Sum Mij nj, with j running over all other levels.
			double creationRate = sourcev(i) + Mvv.row(i) * nv;
			nv(i) = creationRate <= 0 ? 0 : creationRate / fracDestructionRatev(i);
		}

		/* Renormalize because the algorithm has no sum rule, */
		nv *= totalDensity / nv.sum();

		// We will cut this iteration short as long as the ground state has not
		// converged. Use auto to store an expression instead of an intermediate result.
		auto deltaXv = (nv.segment(startX, numX) - previousNv.segment(startX, numX))
		                               .array()
		                               .abs();
		if ((deltaXv > maxDeltaX).any())
		{
			if (!(counter % 1000))
			{
				DEBUG("Solving h2... " << counter << "iterations" << std::endl);
				// DEBUG("h2Levelv = \n" << nv << std::endl);
			}
			continue;
		}

		// If the ground state has more or less converged, we will also start sweeping
		// over the other states. Since theres no dependence between the excited states,
		// we can do this in a block operation (instead of 1 by 1). We use the auto
		// keyword here to work with Eigen expressions, which delays the evaluation
		// without having to put everything on one line.

		// Creation rates for all E levels, coming only from the ground state
		// populations. Calculate intermediate result here, since we will need to check
		// all of them for being > 0 anyway.
		EArray creationRateEv = sourcev.segment(startE, numE) +
		                        Mvv.block(startE, startX, numE, numX) *
		                                        nv.segment(startX, numX);

		// Only the elements for which the denumerator is positive should be calculated,
		// by storing the expression with auto here.
		auto newPopEv = creationRateEv /
		                fracDestructionRatev.segment(startE, numE).array();
		nv.segment(startE, numE) = (creationRateEv > 0).select(newPopEv, 0);

		nv *= totalDensity / nv.sum();

		// Overall convergence check (relative change), also check for changes from 0 to
		// not 0;
		bool fracconv = true;
		for (int i = 0; i < nv.size() && fracconv; i++)
		{
			// finite/0 means not converged
			if (previousNv(i) <= 0 && nv(i) > 0)
				fracconv = false;
			else
			{
				double df = abs(nv(i) / previousNv(i) - 1.);
				if (df > maxFracDelta)
					fracconv = false;
			}
		}

		if (!(counter % 1000))
		{
			DEBUG("Solving h2... " << counter << "iterations" << std::endl);
		}

		converged = fracconv;
	}
	if (nv.array().isNaN().any())
	{
		// TODO: warn about this
		// Error::runtime("nan in H2 level solution");
		nv = (nv.array().isNaN()).select(0, nv);
	}
	DEBUG("Solved H2 in " << counter << " iterations");
	if (converged)
		DEBUG("\n");
	else
		DEBUG(" (not converged)\n");
	// DEBUG("h2Levelv = \n" << nv << std::endl);
	return nv;
}

EVector LevelSolver::statisticalEquilibrium_boltzman(double totalDensity, double T,
                                                     const EVector& ev, const EVector& gv)
{
	// It would be better if we could assume that the lowest level is index 0, but I'm not
	// sure anymore if I made this guaranteed throughout the code.
	double eMin = ev.minCoeff();
	double kT = Constant::BOLTZMAN * T;

	                // g * exp((eMin - e) / kT)
	                EVector pv = gv.array() * ((eMin - ev.array()) / kT).exp();

	return pv / pv.sum() * totalDensity;
}
