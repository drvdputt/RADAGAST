#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"

#include <iostream>

ChemistrySolver::ChemistrySolver(const EMatrix& reactantStoichvv, const EMatrix& productStoichvv,
                                 const EMatrix& conservationCoeffvv)
                : _rvv(reactantStoichvv), _pvv(productStoichvv), _cvv(conservationCoeffvv)
{
	_numSpecies = _rvv.rows();
	_numReactions = _rvv.cols();
	_numConserved = _cvv.rows();
}

ChemistrySolver::ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn) : _cn(std::move(cn))
{
	_rvv = _cn->reactantStoichvv();
	_pvv = _cn->productStoichvv();
	_cvv = _cn->conservationCoeffvv();
	_numSpecies = _rvv.rows();
	_numReactions = _rvv.cols();
	_numConserved = _cvv.rows();
}

ChemistrySolver::~ChemistrySolver() = default;

EVector ChemistrySolver::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
	EVector conservedQuantityv = _cvv * n0v;
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
	fv.head(_numSpecies) = (_pvv - _rvv) * kTotalv;

	//-------------------------------------//
	// BOTTOM PART: CONSERVATION EQUATIONS //
	//-------------------------------------//

	/* These equations are of the form CONSERV_COEFF * DENSITY - CONSTANT = 0. _cvv is indexed
	   on (conserved quantity, species). */
	fv.tail(_numConserved) = _cvv * nv - conservedQuantityv;

	return fv;
}

EMatrix ChemistrySolver::evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const
{
	/* The elements are d f_i / d n_j, so the jacobian should have (dimension of f, dimension of
	   nv). */
	EMatrix jvv(_numSpecies + _numConserved, _numSpecies);

	EMatrix netStoichvv = _pvv - _rvv;

	// For every column (= derivative with respect to a different density)
	for (size_t j = 0; j < _numSpecies; j++)
	{
		/* Here we calculate the derivatives of the density products times the reaction
		   rates and multiply with the rate coefficient. */
		EVector kDerivativev(_numReactions);
		for (size_t r = 0; r < _numReactions; r++)
			kDerivativev(r) = rateCoeffv(r) * densityProductDerivative(nv, r, j);

		// Fill in the top part of the column (= equilibrium part)
		jvv.col(j).head(_numSpecies) = netStoichvv * kDerivativev;

		/* Fill in the bottom part of the column (= conservation part) These equations are
		   simply linear, therefore the derivative is equal to the coefficient. */
		jvv.col(j).tail(_numConserved) = _cvv.col(j);
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
		double stoichRs = _rvv(s, r);
		if (stoichRs)
			densityProduct *= pow(nv(s), stoichRs);
	}
	return densityProduct;
}

double ChemistrySolver::densityProductDerivative(const EVector& nv, size_t r, size_t j) const
{
	/* If the reactant j is not present in reaction (remember that _rvv contains the
	   stoichiometry of the reactants), deriving with respect to it will produce zero. */
	double stoichRj = _rvv(j, r);
	if (!stoichRj)
		return 0;
	else
	{
		// Derivative of the n_j factor
		double densityProductDerivative = stoichRj * pow(nv(j), stoichRj - 1);

		// The rest of the factors
		for (size_t s = 0; s < _numSpecies; s++)
		{
			double stoichRs = _rvv(s, r);
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
	int counter = 0;
	EVector xv = x0v;
	bool converged = true;
	do
	{
		counter++;
		const EMatrix& jfvv = jacobianfvv(xv);
		const EVector& fv = functionv(xv);
		// TODO: remove rows (= equations) that are zero,
		std::cout << "Newton-Raphson iteration " << counter << std::endl;
		std::cout << "Jacobian: " << std::endl << jfvv << std::endl;
		std::cout << "Function: " << std::endl << fv << std::endl;
		EVector deltaxv = jfvv.colPivHouseholderQr().solve(-fv);

		/* Force the result to be positive: x + deltax >= 0 <==> deltax >= -x ==> deltax =
		   max(deltax, -x). */
		deltaxv = deltaxv.array().max(-xv.array());

		// TODO: use 'line search' instead to enforce positivity

		/* We assume that a density has converged when the relative change is less than .1%,
		   or when its relative abundance is negligibly small. Notice that the inequalities
		   NEED to be 'smaller than or equal', in case one of the x'es is zero. */
		Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
		                (deltaxv.array().abs() <= 1.e-17 * xv.array().abs()) ||
		                xv.array() <= 1.e-99 * xv.norm();
		converged = convergedv.all() || (fv.array() == 0).all();

		std::cout << "Delta x:\n"
		          << deltaxv << std::endl
		          << "previous x:\n"
		          << xv << std::endl;
		std::cout << "Convergedv: " << std::endl << convergedv << std::endl;

		xv += deltaxv;
	} while (!converged);
	return xv;
}
