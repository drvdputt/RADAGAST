#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"
#include "DebugMacros.h"
#include "Error.h"

#include <algorithm>
#include <iostream>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sort.h>

ChemistrySolver::ChemistrySolver(const EMatrix& reactantStoichvv,
                                 const EMatrix& productStoichvv,
                                 const EMatrix& conservationCoeffvv)
                : _rStoichvv(reactantStoichvv),
                  _netStoichvv(productStoichvv - reactantStoichvv),
                  _conservEqvv(conservationCoeffvv)
{
	_numSpecies = _rStoichvv.rows();
	_numReactions = _rStoichvv.cols();
	_numConserved = _conservEqvv.rows();
}

ChemistrySolver::ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn) : _cn(std::move(cn))
{
	_rStoichvv = _cn->reactantStoichvv();
	_netStoichvv = _cn->productStoichvv() - _rStoichvv;
	_conservEqvv = _cn->conservationCoeffvv();
	_numSpecies = _rStoichvv.rows();
	_numReactions = _rStoichvv.cols();
	_numConserved = _conservEqvv.rows();
}

ChemistrySolver::~ChemistrySolver() = default;

namespace
{
struct system_params
{
	const EVector* rateCoeffv;
	const EVector* conservedQuantityv;
	const std::vector<size_t>* replaceByConservationv;
	const ChemistrySolver* solver_instance;
};

int system_f(const gsl_vector* x, void* p, gsl_vector* f)
{
	auto* params = static_cast<struct system_params*>(p);

	// View the gsl vector as an eigen vector
	Eigen::Map<EVector> nv(x->data, x->size);

	// Calculate Fv using eigen library
	EVector fv = params->solver_instance->evaluateFv(nv, *params->rateCoeffv,
	                                                 *params->conservedQuantityv,
	                                                 *params->replaceByConservationv);

	// View the resulting eigen vector as a gsl vector, and copy the results
	gsl_vector_view f_view = gsl_vector_view_array(fv.data(), fv.size());
	gsl_vector_memcpy(f, &f_view.vector);
	return GSL_SUCCESS;
}

int system_df(const gsl_vector* x, void* p, gsl_matrix* J)
{
	auto* params = static_cast<system_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);

	// Calculate Jvv using this member function
	EMatrix Jvv = params->solver_instance->evaluateJvv(nv, *params->rateCoeffv,
	                                                   *params->replaceByConservationv);

	// Remember that Eigen works column major, while gsl matrices are indexed row major. So
	// we start by taking a view of the transpose. Should be square though, so the size
	// arguments actually doesn't matter.
	gsl_matrix_view J_view_transposed =
	                gsl_matrix_view_array(Jvv.data(), Jvv.cols(), Jvv.rows());
	gsl_matrix_transpose_memcpy(J, &J_view_transposed.matrix);
	return GSL_SUCCESS;
}

int system_fdf(const gsl_vector* x, void* p, gsl_vector* f, gsl_matrix* J)
{
	auto* params = static_cast<system_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);
	EVector fv = params->solver_instance->evaluateFv(nv, *params->rateCoeffv,
	                                                 *params->conservedQuantityv,
	                                                 *params->replaceByConservationv);
	EMatrix Jvv = params->solver_instance->evaluateJvv(nv, *params->rateCoeffv,
	                                                   *params->replaceByConservationv);
	gsl_vector_view f_view = gsl_vector_view_array(fv.data(), fv.size());
	gsl_vector_memcpy(f, &f_view.vector);
	gsl_matrix_view J_view_transposed =
	                gsl_matrix_view_array(Jvv.data(), Jvv.cols(), Jvv.rows());
	gsl_matrix_transpose_memcpy(J, &J_view_transposed.matrix);
	return GSL_SUCCESS;
}
} /* namespace */

EVector ChemistrySolver::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
#define ROOT_METHOD
#ifdef ROOT_METHOD
	// Fix the conserved quantities using the initial condition
	EVector conservedQuantityv = _conservEqvv * n0v;

	// We will need to replace some of the equations by conservation equations. Numerical
	// stability can be improved by choosing the equations that need to be replaced
	// dynamically. The 2006 meudon paper (section 7) recommends that, for each atom's
	// conservation equation, the equation for the species that is most abundant containing
	// this specific atom is replaced.
	std::vector<size_t> iToReplacev(0);

	// Conserved quantities and rate coefficients are supplied to the functions via this
	// struct
	struct system_params params = {&rateCoeffv, &conservedQuantityv, &iToReplacev, this};

	const gsl_multiroot_fdfsolver_type* T = gsl_multiroot_fdfsolver_hybridsj;
	gsl_multiroot_fdfsolver* s = gsl_multiroot_fdfsolver_alloc(T, _numSpecies);

	gsl_multiroot_function_fdf fdf;
	fdf.f = &system_f;
	fdf.df = &system_df;
	fdf.fdf = &system_fdf;
	fdf.n = _numSpecies;
	fdf.params = &params;

	// Convert the initial guess vector to a gsl_vector without copying
	gsl_vector_const_view initial_guess_view =
	                gsl_vector_const_view_array(n0v.data(), n0v.size());
	gsl_multiroot_fdfsolver_set(s, &fdf, &initial_guess_view.vector);

	double epsabs = 1.e-32;
	double epsrel = 1.e-17;
	// solve here, using maximum 100 iterations
	size_t i = 0;
	gsl_vector* x = gsl_multiroot_fdfsolver_root(s);
	iToReplacev = chooseEquationsToReplace(Eigen::Map<EVector>(x->data, x->size));
	for (; i < 100; i++)
	{
		// Remember that the params struct has a pointer to iToReplacev. The reason
		// (numerical stability) that we need to dynamically choose these indices can be
		// illustraed using a practical example: The H derivative has contributions by
		// H+ recombination and H ionization (fast) and H2 formation and dissociation
		// (slow). If the conservation equations are satisfied (= 0) and we use them to
		// overwrite the H2 equation, then our only knowledge about the formation /
		// dissociation balance is in the H equation. Since these rates are so much
		// smaller than the ones relating to ionization, we effectively lose this
		// information due to precision reasons. By carefully choosing which equations
		// we replace by conservation constraints, we can mitigate problems like this.

		int iterationStatus = gsl_multiroot_fdfsolver_iterate(s);
		if (iterationStatus == GSL_EBADFUNC)
			Error::runtime("Inf or NaN encountered by GSL");
		if (iterationStatus == GSL_ENOPROG)
			iToReplacev = chooseEquationsToReplace(
			                Eigen::Map<EVector>(x->data, x->size));
		// Error::runtime("GSL not making progress");

		x = gsl_multiroot_fdfsolver_root(s);

		gsl_vector* dx = gsl_multiroot_fdfsolver_dx(s);
		gsl_vector* f = gsl_multiroot_fdfsolver_f(s);
		int testDelta = gsl_multiroot_test_delta(dx, x, epsabs, epsrel);
		int testResidual = gsl_multiroot_test_residual(f, epsabs);
		if (testDelta == GSL_SUCCESS && testResidual == GSL_SUCCESS)
		{
			break;
		}
		else if (i == 99)
			Error::runtime("Not converged");
	}

	// Copy the solution
	EVector solutionv(Eigen::Map<EVector>(x->data, x->size));
	auto f = evaluateFv(solutionv, rateCoeffv, conservedQuantityv, std::vector<size_t>{});
	DEBUG(i << " iterations for chemistry\nResidual:\n" << f << '\n');

	gsl_multiroot_fdfsolver_free(s);
#else
	// ode integration or minimization method
#endif
	if (solutionv.array().isNaN().any())
		Error::runtime("NaN in chemistry solution");
	return solutionv;
}

std::vector<size_t> ChemistrySolver::chooseEquationsToReplace(const EVector& nv) const
{
	std::vector<size_t> equationsToReplacev(_numConserved);
	auto iFirst = equationsToReplacev.begin();
	for (int c = 0; c < _numConserved; c++)
	{
		auto iCurrent = iFirst + c;

		// Look at the individual contributions to the conserved quantity (element-wise)
		// Example:

		// species: e- H+ H H2
		// e conservation: 1 0 1 2
		// p conservation: 0 1 1 2
		// contribution = conservation * species density
		EArray contributionv = _conservEqvv.row(c).cwiseProduct(nv.transpose());

		// Do an argsort to sort the species by the contribution size
		std::vector<size_t> speciesIndexv = TemplatedUtils::argsort(
		                &contributionv(0), contributionv.size());

		// Start with the biggest contribution, and go to the next equation (c) when a
		// candidate has been chosen
		for (auto is = speciesIndexv.rbegin(); is != speciesIndexv.rend(); is++)
		{
			// If not already there, choose current candidate index, and move on to
			// the next conservation equation (c)
			if (std::find(iFirst, iCurrent, *is) == iCurrent)
			{
				*iCurrent = *is;
				break;
			}
		}
	}
	return equationsToReplacev;
}

EVector ChemistrySolver::evaluateFv(const EVector& nv, const EVector& rateCoeffv,
                                    const EVector& conservedQuantityv,
                                    const std::vector<size_t>& replaceByConservationv) const
{
	EVector fv(_numSpecies);

	//-------------------------------------------------------------------//
	// TOP PART: EQUILIBRIUM EQUATIONS (A.K.A. NET_STOICH * RATE = ZERO) //
	//-------------------------------------------------------------------//

	// The total rate of each reaction, a.k.a. k(T) multiplied with the correct density
	// factor n_s ^ R_s,r (see notes).
	EVector kTotalv(_numReactions);
	for (int r = 0; r < _numReactions; r++)
		kTotalv(r) = rateCoeffv(r) * densityProduct(nv, r);

	/* The top part of the f-vector is a matrix multiplication between the net stoichiometry
	   matrix (indexed on species, reaction) and the rate vector (indexed on reaction). */
	fv = _netStoichvv * kTotalv;

	//-------------------------------------//
	// BOTTOM PART: CONSERVATION EQUATIONS //
	//-------------------------------------//

	/* These equations are of the form CONSERV_COEFF * DENSITY - CONSTANT = 0. _cvv is
	   indexed on (conserved quantity, species). */
	if (replaceByConservationv.size() > 0)
	{
		Error::equalCheck<int>(
		                "number of conservation equations and list of indices to "
		                "replace",
		                _numConserved, replaceByConservationv.size());
		EVector conservationFv = _conservEqvv * nv - conservedQuantityv;
		// Find the _numConserved smallest rates, and replace their equations by
		// the conservation laws
		for (int i = 0; i < replaceByConservationv.size(); i++)
			fv(replaceByConservationv[i]) = conservationFv(i);
		// fv.tail(_numConserved) = subFv;
	}

	return fv;
}

EMatrix ChemistrySolver::evaluateJvv(const EVector& nv, const EVector& rateCoeffv,
                                     const std::vector<size_t>& replaceByConservationv) const
{
	/* The elements are d f_i / d n_j, so the jacobian should have (dimension of f,
	   dimension of nv). */
	// EMatrix jvv(_numSpecies + _numConserved, _numSpecies);
	// Replace equations
	EMatrix jvv(_numSpecies, _numSpecies);

	// For every column (= derivative with respect to a different density)
	for (size_t jDeriv = 0; jDeriv < _numSpecies; jDeriv++)
	{
		/* Here we calculate the derivatives of the density products times the reaction
		   rates and multiply with the rate coefficient. */
		EVector kDerivativev(_numReactions);
		for (int r = 0; r < _numReactions; r++)
			kDerivativev(r) =
			                rateCoeffv(r) * densityProductDerivative(nv, r, jDeriv);

		jvv.col(jDeriv) = _netStoichvv * kDerivativev;

		// Now overwrite some rows with conservation equations. These equations are
		// simply linear, therefore the derivative is equal to the coefficient.
		// jvv.col(jDeriv).tail(_numConserved) = _conservEqvv.col(jDeriv);
		if (replaceByConservationv.size() > 0)
		{
			Error::equalCheck<int>("number of conservation equations and list of "
			                       "indices to "
			                       "replace",
			                       _numConserved, replaceByConservationv.size());

			// Replace the rows listed in replaceByConservationv
			for (int i = 0; i < replaceByConservationv.size(); i++)
				jvv.col(jDeriv)(replaceByConservationv[i]) =
				                _conservEqvv.col(jDeriv)[i];
			// fv.tail(_numConserved) = subFv;
		}
	}
	return jvv;
}

double ChemistrySolver::densityProduct(const EVector& nv, size_t r) const
{
	double densityProduct = 1;
	for (size_t s = 0; s < _numSpecies; s++)
	{
		/* If the species in involved in this reaction (stoich on left side > 0),
		   calculate the density to the power of its stoichiometry in the reaction. */
		double stoichRs = _rStoichvv(s, r);
		if (stoichRs)
			densityProduct *= pow(nv(s), stoichRs);
	}
	return densityProduct;
}

double ChemistrySolver::densityProductDerivative(const EVector& nv, int r, int j) const
{
	/* If the reactant j is not present in reaction (remember that _rvv contains the
	   stoichiometry of the reactants), deriving with respect to it will produce zero. */
	double stoichRj = _rStoichvv(j, r);
	if (!stoichRj)
		return 0;
	else
	{
		// Derivative of the n_j factor
		double densityProductDerivative = 1;
		if (stoichRj > 1)
			densityProductDerivative *= stoichRj * pow(nv(j), stoichRj - 1);

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
