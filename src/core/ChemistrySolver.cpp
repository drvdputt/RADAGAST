#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"
#include "DebugMacros.h"
#include "Error.h"

#include <algorithm>
#include <iostream>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sort.h>

namespace
{
// Precision
constexpr double epsabs_f = 1e-30;
constexpr double epsabs_x = 1e-30;
constexpr double epsrel_x = 1e-15;

// Stuff for GSL
struct multiroot_params
{
	const EVector* rateCoeffv;
	const EVector* conservedQuantityv;
	const std::vector<size_t>* replaceByConservationv;
	const ChemistrySolver* solver_instance;
};

int multiroot_f(const gsl_vector* x, void* p, gsl_vector* f)
{
	auto* params = static_cast<struct multiroot_params*>(p);

	// View the gsl vector as an eigen vector
	Eigen::Map<EVector> nv(x->data, x->size);

	// if (nv.array().isNaN().any())
	// 	Error::runtime("x is nan");

	// Calculate Fv using eigen library
	EVector fv = params->solver_instance->evaluateFv(nv, *params->rateCoeffv);

	/* Replace some of the Fv by conservation equations, as indicated by the vector of
	   indices. These equations are of the form CONSERV_COEFF * DENSITY - CONSTANT = 0. _cvv
	   is indexed on (conserved quantity, species). */
	const std::vector<size_t>& replaceByConservationv = *params->replaceByConservationv;
	const EVector& q0v = *params->conservedQuantityv;
	if (replaceByConservationv.size() > 0)
	{
		Error::equalCheck<int>(
		                "number of conservation equations and list of indices to "
		                "replace",
		                q0v.size(), replaceByConservationv.size());
		EVector conservationFv = params->solver_instance->evaluateQv(nv, q0v);
		for (int c = 0; c < replaceByConservationv.size(); c++)
			fv(replaceByConservationv[c]) = conservationFv(c);
	}

	// View the resulting eigen vector as a gsl vector, and copy the results
	gsl_vector_view f_view = gsl_vector_view_array(fv.data(), fv.size());
	gsl_vector_memcpy(f, &f_view.vector);
	return GSL_SUCCESS;
}

int multiroot_df(const gsl_vector* x, void* p, gsl_matrix* J)
{
	auto* params = static_cast<multiroot_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);

	// Calculate Jvv using this member function
	EMatrix jvv = params->solver_instance->evaluateJvv(nv, *params->rateCoeffv);

	// Now overwrite some rows with conservation equations. These equations are simply
	// linear, therefore the derivative is equal to the coefficient.
	// jvv.col(jDeriv).tail(_numConserved) = _conservEqvv.col(jDeriv);
	const std::vector<size_t>& replaceByConservationv = *params->replaceByConservationv;
	if (replaceByConservationv.size() > 0)
	{
		Error::equalCheck<int>("number of conservation equations and list of "
		                       "indices to "
		                       "replace",
		                       params->conservedQuantityv->size(),
		                       replaceByConservationv.size());
		EMatrix conservationJvv = params->solver_instance->conservEqvv();
		for (int c = 0; c < replaceByConservationv.size(); c++)
			jvv.row(replaceByConservationv[c]) = conservationJvv.row(c);
	}

	// Remember that Eigen works column major, while gsl matrices are indexed row major. So
	// we start by taking a view of the transpose. Should be square though, so the size
	// arguments actually doesn't matter.
	gsl_matrix_view J_view_transposed =
	                gsl_matrix_view_array(jvv.data(), jvv.cols(), jvv.rows());
	gsl_matrix_transpose_memcpy(J, &J_view_transposed.matrix);
	return GSL_SUCCESS;
}

int multiroot_fdf(const gsl_vector* x, void* p, gsl_vector* f, gsl_matrix* J)
{
	auto* params = static_cast<multiroot_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);

	EVector fv = params->solver_instance->evaluateFv(nv, *params->rateCoeffv);
	EMatrix jvv = params->solver_instance->evaluateJvv(nv, *params->rateCoeffv);

	const std::vector<size_t>& replaceByConservationv = *params->replaceByConservationv;
	if (replaceByConservationv.size() > 0)
	{
		EVector conservationFv = params->solver_instance->evaluateQv(
		                nv, *params->conservedQuantityv);
		EMatrix conservationJvv = params->solver_instance->conservEqvv();
		for (int c = 0; c < replaceByConservationv.size(); c++)
		{
			fv(replaceByConservationv[c]) = conservationFv(c);
			jvv.row(replaceByConservationv[c]) = conservationJvv.row(c);
		}
	}

	gsl_vector_view f_view = gsl_vector_view_array(fv.data(), fv.size());
	gsl_vector_memcpy(f, &f_view.vector);
	gsl_matrix_view J_view_transposed =
	                gsl_matrix_view_array(jvv.data(), jvv.cols(), jvv.rows());
	gsl_matrix_transpose_memcpy(J, &J_view_transposed.matrix);
	return GSL_SUCCESS;
}

struct multimin_params
{
	const EVector* rateCoeffv;
	const EVector* conservedQuantityv;
	const ChemistrySolver* solver_instance;
};

double multimin_f(const gsl_vector* x, void* p)
{
	auto* params = static_cast<struct multimin_params*>(p);

	// View the gsl vector as an Eigen vector (without copying)
	Eigen::Map<EVector> nv(x->data, x->size);

	// Get the vector of derivatives (pass an emtpy vector because we won't be replacing any
	// equations here)
	auto thisp = params->solver_instance;
	EVector Fv = thisp->evaluateFv(nv, *params->rateCoeffv);
	EVector Qv = thisp->evaluateQv(nv, *params->conservedQuantityv);
	double f = thisp->toMinimizeFunction(Fv, Qv);

	// Add penalty for negative densities
	size_t iMin = gsl_vector_min_index(x);
	double min = gsl_vector_get(x, iMin);
	if (min < 0)
		f -= min;
	return f;
}

void multimin_df(const gsl_vector* x, void* p, gsl_vector* g)
{
	auto* params = static_cast<struct multimin_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);

	auto thisp = params->solver_instance;
	EVector Fv = thisp->evaluateFv(nv, *params->rateCoeffv);
	EVector Qv = thisp->evaluateQv(nv, *params->conservedQuantityv);
	EMatrix Jvv = thisp->evaluateJvv(nv, *params->rateCoeffv);

	EVector Gv = thisp->toMinimizeGradientv(Fv, Jvv, Qv);

	// Add slope of penalty for negative densities
	size_t iMin = gsl_vector_min_index(x);
	double min = gsl_vector_get(x, iMin);
	if (min < 0)
		Gv(iMin) -= 1;

	gsl_vector_view g_view = gsl_vector_view_array(Gv.data(), Gv.size());
	gsl_vector_memcpy(g, &g_view.vector);
}

void multimin_fdf(const gsl_vector* x, void* p, double* f, gsl_vector* g)
{
	auto* params = static_cast<struct multimin_params*>(p);
	Eigen::Map<EVector> nv(x->data, x->size);

	auto thisp = params->solver_instance;
	EVector Fv = thisp->evaluateFv(nv, *params->rateCoeffv);
	EVector Qv = thisp->evaluateQv(nv, *params->conservedQuantityv);
	EMatrix Jvv = thisp->evaluateJvv(nv, *params->rateCoeffv);

	*f = thisp->toMinimizeFunction(Fv, Qv);

	EVector Gv = thisp->toMinimizeGradientv(Fv, Jvv, Qv);

	// Add (slope of) penalty for negative densities
	size_t iMin = gsl_vector_min_index(x);
	double min = gsl_vector_get(x, iMin);
	if (min < 0)
	{
		*f -= min;
		Gv(iMin) -= 1;
	}

	gsl_vector_view g_view = gsl_vector_view_array(Gv.data(), Gv.size());
	gsl_vector_memcpy(g, &g_view.vector);
}

struct ode_params
{
	size_t size;
	const EVector* rateCoeffv;
	const ChemistrySolver* solver_instance;
};

int ode_f(double /* unused t */, const double y[], double dydt[], void* p)
{
	auto* params = static_cast<struct ode_params*>(p);

	Eigen::Map<const EVector> nv(y, params->size);
	Eigen::Map<EVector> Fv(dydt, params->size);
	Fv = params->solver_instance->evaluateFv(
	                nv, *params->rateCoeffv); // this should overwrite the contents of dydt
	return GSL_SUCCESS; // maybe check for nan here
}

int ode_j(double /* unused t */, const double y[], double* dfdy, double dfdt[], void* p)
{
	auto* params = static_cast<struct ode_params*>(p);

	// The jacobian dfdy needs to be stored row-major for GSL. This is also why we needed to
	// transpose Jvv in the functions for the multiroot-based solution.
	Eigen::Map<const EVector> nv(y, params->size);
	Eigen::Map<EMatrixRM> Jvv(dfdy, params->size, params->size);
	Jvv = params->solver_instance->evaluateJvv(nv, *params->rateCoeffv);

	// no explicit time dependence
	std::fill(dfdt, dfdt + params->size, 0);

	return GSL_SUCCESS;
}

} /* namespace */

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

EVector ChemistrySolver::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
	// EVector nv = n0v;

	// // First, use the multiroot method, which intentionally uses a slightly modified kv to
	// // prevent the jacobian from becoming singular (resulting in weird steps)
	// double smallest = rateCoeffv.maxCoeff();
	// for (int i = 0; i < rateCoeffv.size(); i++)
	// 	if (rateCoeffv(i) > 0)
	// 		smallest = std::min(smallest, rateCoeffv(i));

	// // Replace the zeros by something very small
	// EVector kv = (rateCoeffv.array() > 0).select(rateCoeffv, smallest);
	// nv = solveMultiroot(kv, nv);

	// // Then, use the minimization algorithm which does use the correct coefficients to go to
	// // the real solution
	// nv = solveMultimin(rateCoeffv, nv);

	// if (nv.array().isNaN().any())
	// 	Error::runtime("NaN in chemistry solution");

	// // Remove negative densities (and hope everything is still normalized)
	// return (nv.array() < 0).select(0, nv);

	return solveTimeDep(rateCoeffv, n0v);
}

EVector ChemistrySolver::solveMultiroot(const EVector& rateCoeffv, const EVector& n0v) const
{
	// Fix the conserved quantities using the initial condition
	EVector conservedQuantityv = _conservEqvv * n0v;

	// We will need to replace some of the equations by conservation equations. Numerical
	// stability can be improved by choosing the equations that need to be replaced
	// dynamically. The 2006 meudon paper (section 7) recommends that, for each atom's
	// conservation equation, the equation for the species that is most abundant containing
	// this specific atom is replaced.
	std::vector<size_t> iToReplacev = chooseEquationsToReplace(n0v);

	// Conserved quantities and rate coefficients are supplied to the functions via this
	// struct
	struct multiroot_params params = {&rateCoeffv, &conservedQuantityv, &iToReplacev, this};

	const gsl_multiroot_fdfsolver_type* T = gsl_multiroot_fdfsolver_hybridsj;
	gsl_multiroot_fdfsolver* s = gsl_multiroot_fdfsolver_alloc(T, _numSpecies);

	gsl_multiroot_function_fdf fdf;
	fdf.f = &multiroot_f;
	fdf.df = &multiroot_df;
	fdf.fdf = &multiroot_fdf;
	fdf.n = _numSpecies;
	fdf.params = &params;

	// Convert the initial guess vector to a gsl_vector without copying
	gsl_vector_const_view initial_guess_view =
	                gsl_vector_const_view_array(n0v.data(), n0v.size());
	gsl_multiroot_fdfsolver_set(s, &fdf, &initial_guess_view.vector);

	// solve here, using maximum 100 iterations
	size_t i = 0;
	gsl_vector* x = gsl_multiroot_fdfsolver_root(s);
	gsl_vector* dx = gsl_multiroot_fdfsolver_dx(s);
	gsl_vector* f = gsl_multiroot_fdfsolver_f(s);
	int maxIt = 100;
	for (; i < maxIt; i++)
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
		iToReplacev = chooseEquationsToReplace(Eigen::Map<EVector>(x->data, x->size));

		int iterationStatus = gsl_multiroot_fdfsolver_iterate(s);
		int testDelta = gsl_multiroot_test_delta(dx, x, epsabs_x, epsrel_x);
		int testResidual = gsl_multiroot_test_residual(f, epsabs_f);

		if (testDelta == GSL_SUCCESS && testResidual == GSL_SUCCESS)
			break;
		if (iterationStatus == GSL_ENOPROG)
		{
			DEBUG("No progress! breaking...\n");
			break;
		}
		if (iterationStatus == GSL_EBADFUNC)
			Error::runtime("Inf or NaN encountered by GSL");
	}

	// Copy the solution
	EVector solutionv(Eigen::Map<EVector>(x->data, x->size));
	EVector Fv = evaluateFv(solutionv, rateCoeffv);
	EVector Qv = evaluateQv(solutionv, conservedQuantityv);
	DEBUG(i << "multiroot iterations for chemistry\nResidual:\n"
	        << Fv << "\nConservation\n"
	        << Qv << "\nSolution\n"
	        << solutionv << "\n");

	gsl_multiroot_fdfsolver_free(s);
	return solutionv;
}

EVector ChemistrySolver::solveMultimin(const EVector& rateCoeffv, const EVector& n0v) const
{
	EVector conservedQuantityv = _conservEqvv * n0v;

	struct multimin_params params = {&rateCoeffv, &conservedQuantityv, this};

	const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_conjugate_fr;
	gsl_multimin_fdfminimizer* m = gsl_multimin_fdfminimizer_alloc(T, _numSpecies);

	gsl_multimin_function_fdf fdf;
	fdf.f = &multimin_f;
	fdf.df = &multimin_df;
	fdf.fdf = &multimin_fdf;
	fdf.n = _numSpecies;
	fdf.params = &params;

	gsl_vector_const_view initial_guess_view =
	                gsl_vector_const_view_array(n0v.data(), n0v.size());
	double step_size = n0v.norm();
	double tol = 0.1;
	gsl_multimin_fdfminimizer_set(m, &fdf, &initial_guess_view.vector, step_size, tol);

	double epsabs_f2 = epsabs_f * epsabs_f;
	double epsabs_g = epsabs_f2 / epsabs_x;
	int maxIt = 100;
	gsl_vector* x = gsl_multimin_fdfminimizer_x(m);
	double minimum = gsl_multimin_fdfminimizer_minimum(m);
	gsl_vector* g = gsl_multimin_fdfminimizer_gradient(m);

	// Map for easy nan checking
	Eigen::Map<EVector> solutionv(x->data, x->size);

	// Store a fallback solution that was not NaN, in case we get a NaN.
	EVector nonNanSolutionv(solutionv);

	int i = 0;
	for (; i < maxIt; i++)
	{
		int iterationStatus = gsl_multimin_fdfminimizer_iterate(m);
		minimum = gsl_multimin_fdfminimizer_minimum(m);
		bool testMinimum = minimum < epsabs_f2;
		int testGradient = gsl_multimin_test_gradient(g, epsabs_g);

		// DEBUG("Step\n" << Eigen::Map<EVector>(dx->data, dx->size) << '\n');

		if (solutionv.array().isNaN().any())
			break;
		else
			nonNanSolutionv = solutionv;
		if (iterationStatus == GSL_ENOPROG)
		{
			DEBUG("No progress! breaking...\n");
			break;
		}
		if (testGradient == GSL_SUCCESS && testMinimum)
			break;
	}

	EVector Fv = evaluateFv(nonNanSolutionv, rateCoeffv);
	EVector Qv = evaluateQv(nonNanSolutionv, conservedQuantityv);
	DEBUG(i << "multimin iterations for chemistry\nResidual:\n"
	        << Fv << "\nConservation\n"
	        << Qv << "\nSolution\n"
	        << nonNanSolutionv << "\n");

	gsl_multimin_fdfminimizer_free(m);
	return nonNanSolutionv;
}

EVector ChemistrySolver::solveTimeDep(const EVector& rateCoeffv, const EVector& n0v) const
{
	EVector conservedQuantityv = _conservEqvv * n0v;

	struct ode_params params = {_numSpecies, &rateCoeffv, this};

	gsl_odeiv2_system s{ode_f, ode_j, _numSpecies, &params};

	// 1 second initial step
	double ini_step = 1;
	double epsabs = 1;
	double epsrel = 1e-15;
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&s, gsl_odeiv2_step_bsimp,
	                                                     ini_step, epsabs, epsrel);

	// 2^59 seconds > 13.7 Gyr
	int max_steps = 60;
	double t = 0;
	EVector nv = n0v;
	for (int i = 1; i <= max_steps; i++)
	{
		EVector previousNv = nv;

		double goalTime = std::pow(2., i);
		int status = gsl_odeiv2_driver_apply(d, &t, goalTime, &nv[0]);
		if (status != GSL_SUCCESS)
		{
			DEBUG("GSL failure" << '\n');
			// If there was a problem, then nv might have become NaN. Use the last
			// good value.
			nv = previousNv;
			break;
		}
		EArray delta = (nv - previousNv).array().abs();
		EArray avg = (nv + previousNv) / 2;
		bool absEquil = (delta < epsabs).all();
		bool relEquil = ((delta / avg).abs() < epsrel).all();
		if (absEquil && relEquil)
		{
			DEBUG("Reached chemical equilibrium after " << i << " iterations\n");
			break;
		}
	}
	gsl_odeiv2_driver_free(d);
	return nv;
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

EVector ChemistrySolver::evaluateFv(const EVector& nv, const EVector& rateCoeffv) const

{
	EVector fv(_numSpecies);

	//---------------------------------------------------------//
	// EQUILIBRIUM EQUATIONS (A.K.A. NET_STOICH * RATE = ZERO) //
	//---------------------------------------------------------//

	// The total rate of each reaction, a.k.a. k(T) multiplied with the correct density
	// factor n_s ^ R_s,r (see notes).
	EVector kTotalv(_numReactions);
	for (int r = 0; r < _numReactions; r++)
		kTotalv(r) = rateCoeffv(r) * densityProduct(nv, r);

	/* The top part of the f-vector is a matrix multiplication between the net stoichiometry
	   matrix (indexed on species, reaction) and the rate vector (indexed on reaction). */
	fv = _netStoichvv * kTotalv;

	return fv;
}

EMatrix ChemistrySolver::evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const

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
	}
	return jvv;
}

EVector ChemistrySolver::evaluateQv(const EVector& nv, const EVector& conservedQuantityv) const
{
	return _conservEqvv * nv - conservedQuantityv;
}

double ChemistrySolver::toMinimizeFunction(const EVector& Fv, const EVector& Qv) const
{
	return Fv.squaredNorm() + Qv.squaredNorm();
}

EVector ChemistrySolver::toMinimizeGradientv(const EVector& Fv, const EMatrix& Jvv,
                                             const EVector& Qv) const
{
	// I calculated this on paper: gradient of (v.v) is Jv^T v. The jacobian of the
	// conservation equations is just the coefficient matrix _conservEqvv.

	// sizes: //// n x n //////// n //////////// n x c ////////// c //
	// n = numSpecies; c = num of conserved quantities
	return 2 * (Jvv.transpose() * Fv + _conservEqvv.transpose() * Qv);
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
