#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"
#include "DebugMacros.h"

#include <algorithm>
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

EVector ChemistrySolver::newtonRaphson(std::function<EMatrix(const EVector& xv)> functionJvv,
                                       std::function<EVector(const EVector& xv)> functionFv,
                                       const EVector& x0v) const
{
	int iNRIteration{0};
	EVector xv = x0v;
	bool converged = true;
	do
	{
		iNRIteration++;

		// Do a first evaluation of jvv and fv to test the waters
		EMatrix jvv = functionJvv(xv);
		EVector fv = functionFv(xv);

		/* I want to remove species (== 1 row and 1 column) that cause numerical instability
		   A typical culprit is a very high or very low dissociation rate for H2
		   Then, we get like [small small small big] for H and [0 0 0 -big] for the H2 row. */
		std::vector<bool> allZerov(_numSpecies, false);
		std::vector<bool> allZeroButDiagonalv(_numSpecies, false);
		std::vector<bool> isRemovedv(_numSpecies, false);
		int countToRemove{0};
		for (size_t i = 0; i < _numSpecies; i++)
		{
			if ((jvv.row(i).array() == 0).all())
				allZerov[i] = true;
			else
			{
				bool restZero{true};
				for (size_t j = 0; j < _numSpecies && restZero; j++)
					if (j != i && jvv(i, j))
						restZero = false;
				allZeroButDiagonalv[i] == restZero&& jvv(i, i);
			}
			if (allZerov[i] || allZeroButDiagonalv[i])
			{
				isRemovedv[i] = true;
				countToRemove++;
			}
		}

		// Split up the rest of the implementation
		EVector deltaxv(_numSpecies);
		if (!countToRemove)
		{
			// If no equations/species need to be removed, simply take a Newton-Raphson step
			deltaxv = newtonRaphsonStep(jvv, fv, functionFv, xv);
		}
		else
		{
			// If there are some species/equations giving the algorithm trouble, we will redefine the system here
			// No doubt that this is very hacky, and possibly a bit too specialized

			size_t newNumEq = jvv.rows() - countToRemove;
			size_t newNumSp = jvv.cols() - countToRemove;

			/* If all the elements of a row are zero, the density should stay constant
			   If all the elements except the 'diagonal' element are zero, then the density
			   of a species should become zero. */

			/* Function that converts the 'reduced' density vector to a modified 'full' one. This algorithm takes the current value of the full density vector,
			   sets the elements that need to be zero to zero, and puts the values of the reduced density vector in the right place. */
			auto fullXvWithConstants = [&](const EVector& reducedXv) -> EVector {
				EVector result{xv};
				int reducedIndex{0};
				for (size_t i = 0; i < _numSpecies; i++)
				{
					// All coefficients are zero -> keep current value constant
					// *do nothing*
					// Only the diagonal element is nonzero, so density has to be zero
					if (allZeroButDiagonalv[i])
						result(i) = 0;
					// A normal row of coefficients, so the density may evolve
					else if (!allZerov[i])
					{
						result(i) = reducedXv(i);
						reducedIndex++;
					}
				}
				return result;
			};

			/* Define a new jacobian and residual. These functions take a reduced density vector, and internally convert it to a full one. Only the
			   elements that are in the reduced on will evolve when a newton raphson step is taken using these two functions. */
			auto functionReducedJvv = [&](const EVector& reducedXv) -> EMatrix {

				// Evaluate the full jacobian for this modified density vector
				EMatrix wholeJfvv = functionJvv(fullXvWithConstants(reducedXv));

				// Use only the rows and columns which correspond to non-constant species, and the (reduced) conservation equations
				EMatrix result(newNumEq, newNumSp);
				int row{0};
				for (size_t i = 0; i < wholeJfvv.rows(); i++)
					if (i > _numSpecies || !isRemovedv[i])
					{
						int col{0};
						for (size_t j = 0; j < wholeJfvv.cols(); j++)
							if (!(allZerov[j] ||
							      allZeroButDiagonalv[j]))
							{
								result(row, col) = wholeJfvv(i, j);
								col++;
							}
						row++;
					}
				return result;
			};
			auto functionReducedFv = [&](const EVector& reducedXv) -> EVector {

				// Evaluate the full residual
				EVector wholefv = functionFv(fullXvWithConstants(reducedXv));

				// Use only the elements which correspond to non-constant species, and the conservation equations
				EVector result(newNumEq);
				int col{0};
				for (size_t i = 0; i < wholefv.size(); i++)
					if (i > _numSpecies || !isRemovedv[i])
					{
						result(col) = wholefv(i);
						col++;
					}
				return result;
			};

			// Copy over the current xv, except for the species that are no longer part of the system
			EVector currentReducedXv(newNumSp);
			int reducedIndex{0};
			for (size_t i = 0; i < _numSpecies; i++)
				if (!isRemovedv[i])
				{
					currentReducedXv(reducedIndex) = xv(i);
					reducedIndex++;
				}

			// Do a newton raphson step for the reduced density vector
			jvv = functionReducedJvv(currentReducedXv);
			fv = functionReducedFv(currentReducedXv);
			EVector reducedDeltaxv = newtonRaphsonStep(jvv, fv, functionReducedFv,
			                                           currentReducedXv);

			// Fill in the full deltaxv (zeros for species not included)
			reducedIndex = 0;
			for (size_t i = 0; i < _numSpecies; i++)
				if (isRemovedv[i])
					deltaxv(i) = 0;
				else
				{
					deltaxv(i) = reducedDeltaxv(reducedIndex);
					reducedIndex++;
				}
		}

		DEBUG("Newton-Raphson iteration " << iNRIteration << std::endl);
#ifdef PRINT_MATRICES
		DEBUG("Jacobian: " << std::endl << jvv << std::endl);
		DEBUG("Function: " << std::endl << fv << std::endl);
#endif

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

EVector ChemistrySolver::newtonRaphsonStep(const EMatrix& currentJvv, const EMatrix& currentFv,
                                           std::function<EVector(const EVector& nv)> functionv,
                                           const EVector& currentXv) const
{
	EVector deltaxv = currentJvv.colPivHouseholderQr().solve(-currentFv);
	EVector testProduct = -currentJvv * deltaxv;
	DEBUG("testfv \n" << testProduct << std::endl);

	/* A possible problem is that, when the density of a species i is zero, the > 0
	   density criterion prevents steps from being taken whenever delta x_i is negative.

	   Attempt to fix 1: whenever this happens, just set the component delta x_i to
	   zero, and do the line search along the x_i = 0 hyperplane.

	   Result 1: For small dissociation rate (1e-15), and starting with everything H2,
	   a lot of iterations are needed, but
	   the algorithm eventually pushes through. For even smaller rates, numerical
	   instability occurs. */

	/* Never take negative step in x_i when the x_i is zero. */
	for (size_t i = 0; i < currentXv.size(); i++)
		if (currentXv(i) <= 0 && deltaxv(i) < 0)
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
	EVector newxv{currentXv + factor * deltaxv};
	while (count < maxCount)
	{

		/* If all densities are positive, and the norm has decreased, then we have
		   a good step. */
		if ((newxv.array() >= 0).all())
		{
			EVector newfv = functionv(newxv);
			if (newfv.norm() < currentFv.norm())
				break;
		}

		/* In case there's a negative density or an increase of the norm of fv,
		   adjust the step scale factor. */
		factor *= factorReduce;
		newxv = currentXv + factor * deltaxv;
		count++;
	}

	/* If there is still a negative density (count ran out for example), choose the
	   scale factor such that this density is incremented to zero. */
	if ((newxv.array() < 0).any())
	{
		int iMin;
		newxv.array().minCoeff(&iMin);
		// demand that xv(iMin) + factor * deltaxv(iMin) = 0 instead of negative
		factor = -currentXv(iMin) / deltaxv(iMin);
	}

	// Our final choice for deltaxv
	deltaxv *= factor;
	return currentXv + deltaxv;
}
