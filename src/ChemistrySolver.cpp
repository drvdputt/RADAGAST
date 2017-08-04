#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"
#include "DebugMacros.h"

#include <iostream>

ChemistrySolver::ChemistrySolver(const EMatrix& reactantStoichvv, const EMatrix& productStoichvv,
                                 const EMatrix& conservationCoeffvv)
                : _rStoichvv(reactantStoichvv), _netStoichvv(productStoichvv - reactantStoichvv),
                  _conservEqvv(conservationCoeffvv)
{
	_numSpecies = _rStoichvv.rows();
	_numReactions = _rStoichvv.cols();
	_numConserved = _conservEqvv.rows();
}

ChemistrySolver::ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn) : _cn(std::move(cn))
{
	_rStoichvv = _cn->reactantStoichvv();
	_netStoichvv = _rStoichvv - _cn->productStoichvv();
	_conservEqvv = _cn->conservationCoeffvv();
	_numSpecies = _rStoichvv.rows();
	_numReactions = _rStoichvv.cols();
	_numConserved = _conservEqvv.rows();
}

ChemistrySolver::~ChemistrySolver() = default;

EVector ChemistrySolver::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
	EVector conservedQuantityv = _conservEqvv * n0v;
	EVector result = newtonRaphson(
	                [&](const EVector& nv) { return evaluateJvv(nv, rateCoeffv); },
	                [&](const EVector& nv) {
		                return evaluateFv(nv, rateCoeffv, conservedQuantityv);
		        },
	                n0v);
	return result;
}

EVector ChemistrySolver::evaluateFv(const EVector& nv, const EVector& rateCoeffv,
                                    const EVector& conservedQuantityv) const
{
	EVector fv(_numSpecies + _numConserved);

	//-------------------------------------------------------------------//
	// TOP PART: EQUILIBRIUM EQUATIONS (A.K.A. NET_STOICH * RATE = ZERO) //
	//-------------------------------------------------------------------//

	// The total reaction rates, a.k.a. k(T) multiplied with the correct density factors
	EVector kTotalv(_numReactions);

	/* For every reaction, calculate the product of the densities: n_s ^ R_s,r and multiply with
	   the rate coefficient. */
	for (size_t r = 0; r < _numReactions; r++)
		kTotalv(r) = rateCoeffv(r) * densityProduct(nv, r);

	/* The top part of the f-vector is a matrix multiplication between the net stoichiometry
	   matrix (indexed on species, reaction) and the rate vector (indexed on reaction). */
	fv.head(_numSpecies) = _netStoichvv * kTotalv;

	//-------------------------------------//
	// BOTTOM PART: CONSERVATION EQUATIONS //
	//-------------------------------------//

	/* These equations are of the form CONSERV_COEFF * DENSITY - CONSTANT = 0. _cvv is indexed
	   on (conserved quantity, species). */
	fv.tail(_numConserved) = _conservEqvv * nv - conservedQuantityv;

	return fv;
}

EMatrix ChemistrySolver::evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const
{
	/* The elements are d f_i / d n_j, so the jacobian should have (dimension of f, dimension of
	   nv). */
	EMatrix jvv(_numSpecies + _numConserved, _numSpecies);

	// For every column (= derivative with respect to a different density)
	for (size_t j = 0; j < _numSpecies; j++)
	{
		/* Here we calculate the derivatives of the density products times the reaction
		   rates and multiply with the rate coefficient. */
		EVector kDerivativev(_numReactions);
		for (size_t r = 0; r < _numReactions; r++)
			kDerivativev(r) = rateCoeffv(r) * densityProductDerivative(nv, r, j);

		// Fill in the top part of the column (= equilibrium part)
		jvv.col(j).head(_numSpecies) = _netStoichvv * kDerivativev;

		/* Fill in the bottom part of the column (= conservation part) These equations are
		   simply linear, therefore the derivative is equal to the coefficient. */
		jvv.col(j).tail(_numConserved) = _conservEqvv.col(j);
	}
	return jvv;
}

double ChemistrySolver::densityProduct(const EVector& nv, size_t r) const
{
	double densityProduct = 1;
	for (size_t s = 0; s < _numSpecies; s++)
	{
		/* If the species in involved in this reaction (stoich on left side > 0), calculate
		   the density to the power of its stoichiometry in the reaction. */
		double stoichRs = _rStoichvv(s, r);
		if (stoichRs)
			densityProduct *= pow(nv(s), stoichRs);
	}
	return densityProduct;
}

double ChemistrySolver::densityProductDerivative(const EVector& nv, size_t r, size_t j) const
{
	/* If the reactant j is not present in reaction (remember that _rvv contains the
	   stoichiometry of the reactants), deriving with respect to it will produce zero. */
	double stoichRj = _rStoichvv(j, r);
	if (!stoichRj)
		return 0;
	else
	{
		// Derivative of the n_j factor
		double densityProductDerivative = stoichRj * pow(nv(j), stoichRj - 1);

		// The rest of the factors
		for (size_t s = 0; s < _numSpecies; s++)
		{
			double stoichRs = _rStoichvv(s, r);
			if (stoichRs && s != j)
				densityProductDerivative *= pow(nv(s), stoichRs);
		}
		return densityProductDerivative;
	}
}

EVector ChemistrySolver::newtonRaphson(std::function<EMatrix(const EVector& xv)> jacobianfvv,
                                       std::function<EVector(const EVector& xv)> functionv,
                                       const EVector& x0v) const
{
	int iNRIteration{0};
	EVector xv = x0v;
	bool converged = true;
	do
	{
		iNRIteration++;
		const EMatrix& jfvv = jacobianfvv(xv);
		const EVector& fv = functionv(xv);
// NOT YET OPERATIONAL -VVV-
// I ALSO FEEL LIKE THIS IS BECOMING TOO SPECIALIZED. I SHOULD PROBABLY WAIT UNTIL LATER FOR
//		STABILIZING / OPTIMIZING
//		std::vector<size_t> allZerov;
//		std::vector<size_t> allZeroButDiagonalv;
//		for (size_t i = 0; i < _numSpecies; i++)
//		{
//			if ((jfvv.row(i).array() == 0).all())
//				allZerov.emplace_back(i);
//
//			bool restZero{true};
//			for (size_t j = 0; j < _numSpecies && restZero; j++)
//				if (j != i && jfvv(i, j))
//					restZero = false;
//			if (restZero && jfvv(i, i))
//				allZeroButDiagonalv.emplace_back(i);
//		}
//		/* If all the elements of a row are zero, the density should stay constant */
//		/* If all the elements except the 'diagonal' element are zero, then the density
//		   of a species should become zero */
//		size_t numAllZero = allZerov.size();
//		size_t numAllZeroButDiagonal = allZeroButDiagonalv.size();
//		size_t numRemoved = numAllZero + numAllZeroButDiagonal;
//		if (numRemoved)
//		{
//			size_t newNumEq = jfvv.rows() - numRemoved;
//			size_t newNumSp = jfvv.cols() - numRemoved;
//
//			// Define a new jacobian and residual
//			auto newJacobianvv = [&](const EVector& xv) -> EMatrix {
//				EMatrix result(newNumEq, newNumSp);
//				EVector wholeXv(x0v.size());
//
//				EMatrix wholeJfvv = jacobianfvv(xv);
//				return result;
//			};
//			auto newFuntionv = [&](const EVector& xv) > -> EVector {
//				EVector result(newNumEq);
//				return result;
//			};
//
//			// Remove the relevant rows and columns, and keep the species constant
//		}

		/* TODO: remove species (== 1 row and 1 column)  that cause numerical instability
		   A typical culprit is a very high or very low dissociation rate for H2
		   Then, we get like [small small small big] for H and [0 0 0 -big] for the H2 row.

		   Another problem is that, when the density of a species i is zero, the > 0
		   density criterion prevents steps from being taken whenever delta x_i is negative.

		   Attempt to fix 1: whenever this happens, just set the component delta x_i to
		   zero, and do the line search along the x_i = 0 hyperplane.

		   Result 1: For small dissociation rate (1e-15), and starting with everything H2,
		   a lot of iterations are needed, but
		   the algorithm eventually pushes through. For even smaller rates, numerical
		   instability occurs. */

		DEBUG("Newton-Raphson iteration " << iNRIteration << std::endl);
#ifdef PRINT_MATRICES
		DEBUG("Jacobian: " << std::endl << jfvv << std::endl);
		DEBUG("Function: " << std::endl << fv << std::endl);
#endif

		EVector deltaxv = jfvv.colPivHouseholderQr().solve(-fv);
		EVector testProduct = -jfvv * deltaxv;
		DEBUG("testfv \n" << testProduct << std::endl);

		/* Never take negative step in x_i when the x_i is zero. */
		for (size_t i = 0; i < _numSpecies; i++)
			if (xv(i) <= 0 && deltaxv(i) < 0)
				deltaxv(i) = 0;

		/* Line search such as described in Numerical Recipes (I'm basing this on a
		   flowchart found in Comput. Chem. Eng., 2013, 58, 135 - 143). The bottom line is,
		   we rescale the step until none of the densities are negative, and the norm of
		   the residual function ||f(x + factor * deltax)|| will actually become smaller. Note
		   that this will not always be that case for factor=1, as the system is non-linear.
		   The second criterion can
		   be dropped if no solution with a smaller ||f|| can be found. Another option is
		   finding the minimum along the line. I'll keep that in mind for later maybe. */
		double factor{1.}, factorReduce{0.9};
		int count{0}, maxCount{100};
		EVector newxv{xv + factor * deltaxv};
		while (count < maxCount)
		{

			/* If all densities are positive, and the norm has decreased, then we have
			   a good step. */
			if ((newxv.array() >= 0).all())
			{
				EVector newfv = functionv(newxv);
				if (newfv.norm() < fv.norm())
					break;
			}

			/* In case there's a negative density or an increase of the norm of fv,
			   adjust the step scale factor. */
			factor *= factorReduce;
			newxv = xv + factor * deltaxv;
			count++;
		}

		/* If there is still a negative density (count ran out for example), choose the
		   scale factor such that this density is incremented to zero. */
		if ((newxv.array() < 0).any())
		{
			int iMin;
			newxv.array().minCoeff(&iMin);
			// demand that xv(iMin) + factor * deltaxv(iMin) = 0 instead of negative
			factor = -xv(iMin) / deltaxv(iMin);
		}

		// Our final choice for deltaxv
		deltaxv *= factor;

		/* We assume that a density has converged when the relative change is less than .1%,
		   or when its relative abundance is negligibly small. Notice that the inequalities
		   NEED to be 'smaller than or equal', in case one of the x'es is zero. */
		Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
		                deltaxv.array().abs() <= 1.e-9 * xv.array().abs() ||
		                xv.array().abs() <= 1.e-99 * xv.norm();
		converged = convergedv.all() || (fv.array() == 0).all();

#ifdef PRINT_MATRICES
		DEBUG("Delta x:\n"
		      << deltaxv << std::endl
		      << "previous x:\n"
		      << xv << std::endl
		      << "Convergedv: " << std::endl
		      << convergedv << std::endl);
#endif

		// Final choice for updated xv
		xv += deltaxv;

	} while (!converged);

	return xv;
}
