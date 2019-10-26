#include "Chemistry.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "SpeciesIndex.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace
{
struct ode_params
{
	size_t size;
	const EVector* rateCoeffv;
	const Chemistry* chemistry;
};

int ode_f(double /* unused t */, const double y[], double dydt[], void* p)
{
	auto* params = static_cast<struct ode_params*>(p);

	Eigen::Map<const EVector> nv(y, params->size);
	Eigen::Map<EVector> Fv(dydt, params->size);
	Fv = params->chemistry->evaluateFv(
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
	Jvv = params->chemistry->evaluateJvv(nv, *params->rateCoeffv);

	// no explicit time dependence
	std::fill(dfdt, dfdt + params->size, 0);

	return GSL_SUCCESS;
}
} // namespace

void Chemistry::addReaction(const std::string& reactionName,
                            const std::vector<std::string>& reactantNamev,
                            const Array& reactantStoichv,
                            const std::vector<std::string>& productNamev,
                            const Array& productStoichv)
{
	// Give the reaction a number, and put its name in the map
	_reactionIndexm.emplace(reactionName, _reactionIndexm.size());
	_reactionv.emplace_back(
	                SpeciesIndex::makeFullCoefficientv(reactantNamev, reactantStoichv),
	                SpeciesIndex::makeFullCoefficientv(productNamev, productStoichv));
}

void Chemistry::prepareCoefficients()
{
	// TODO: numSpecies needs to be more flexible
	_numSpecies = SpeciesIndex::size();
	_rStoichvv = makeReactantStoichvv();
	_netStoichvv = makeProductStoichvv() - _rStoichvv;
	_numReactions = _rStoichvv.cols();
}

int Chemistry::reactionIndex(const std::string& reactionName) const
{
	return _reactionIndexm.at(reactionName);
}

EVector Chemistry::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
	Error::equalCheck<int>("chemistry coefficients", _numReactions, _reactionv.size());
	Error::equalCheck<int>("chemistry coefficients", _numSpecies, _netStoichvv.rows());
	return solveTimeDep(rateCoeffv, n0v);
}

EVector Chemistry::evaluateFv(const EVector& nv, const EVector& rateCoeffv) const

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

	// Matrix multiplication between the net stoichiometry matrix (indexed on species,
	// reaction) and the rate vector (indexed on reaction).
	fv = _netStoichvv * kTotalv;

	return fv;
}

EMatrix Chemistry::evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const

{
	// The elements are d f_i / d n_j, so the jacobian should have (dimension of fv,
	// dimension of nv).
	EMatrix jvv(_numSpecies, _numSpecies);

	// For every column (= derivative with respect to a different density)
	for (size_t jDeriv = 0; jDeriv < _numSpecies; jDeriv++)
	{
		// Here we calculate the derivatives of the density products times the reaction
		// rates and multiply with the net coefficient.
		EVector kDerivativev(_numReactions);
		for (int r = 0; r < _numReactions; r++)
			kDerivativev(r) =
			                rateCoeffv(r) * densityProductDerivative(nv, r, jDeriv);
		jvv.col(jDeriv) = _netStoichvv * kDerivativev;
	}
	return jvv;
}

EMatrix Chemistry::makeReactantStoichvv() const
{
	EMatrix r(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
	{
		// Each column represents a reaction
		r.col(j) = _reactionv[j]._rv;
	}
	return r;
}

EMatrix Chemistry::makeProductStoichvv() const
{
	EMatrix p(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
	{
		// Each column represents a reaction
		p.col(j) = _reactionv[j]._pv;
	}
	return p;
}

EVector Chemistry::solveTimeDep(const EVector& rateCoeffv, const EVector& n0v) const
{
	struct ode_params params = {_numSpecies, &rateCoeffv, this};

	gsl_odeiv2_system s{ode_f, ode_j, _numSpecies, &params};

	// 1 second initial step
	double ini_step = 1;
	double epsabs = 1;
	double epsrel = 1e-15;
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&s, gsl_odeiv2_step_bsimp,
	                                                     ini_step, epsabs, epsrel);

	int max_steps = 32;
	double t = 0;
	EVector nv = n0v;
	for (int i = 1; i <= max_steps; i++)
	{
		// We need to advance time by enough to see a significant difference in the
		// slowest changing densities. If the chosen time scale is too short, then we
		// would wrongly conclude that some densities have reached equilibrium, while in
		// fact they just wouldn't have had enough time to change. Time scale =
		// max(density / rate of change). Don't forget abs because rate of change can be
		// negative of course.

		EArray fv = evaluateFv(nv, rateCoeffv);
		double timeScale = (fv != 0).select(nv.array() / fv.abs(), 0).maxCoeff();

		// Limit the time scale to roughly the age of the universe
		timeScale = std::min(timeScale, 5.5e17);

		// I found that using twice the value works slightly better, but there's a lot
		// of wiggle room of course
		double goalTime = t + 2 * timeScale;

		EVector previousNv = nv;
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
		// Ignore relative change when denominator is zero
		bool relEquil = ((delta / avg).abs() < epsrel || avg == 0).all();
		if (absEquil && relEquil)
		{
			DEBUG("Reached chemical equilibrium after " << i << " iterations\n");
			break;
		}
	}
	gsl_odeiv2_driver_free(d);
	return nv;
}

double Chemistry::densityProduct(const EVector& nv, size_t r) const
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

// TODO: speed this up by calculating the whole jacobian of kTotalv at once. (Powers nj^Rjr
// often recurr between different j). Should also be easily testable).
double Chemistry::densityProductDerivative(const EVector& nv, int r, int j) const
{
	// If the reactant j is not present in reaction (remember that _rvv contains the
	// stoichiometry of the reactants), deriving with respect to it will produce zero. */
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
