#include "Chemistry.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

namespace
{
struct ode_params
{
	size_t size;
	const EVector* rateCoeffv;
	const Chemistry* chemistry;
	// Workspace for total rate vector (prevent frequent reallocation)
	EVector* kv;
	EMatrix* Jkvv;
};

int ode_f(double /* unused t */, const double y[], double dydt[], void* p)
{
	auto* params = static_cast<struct ode_params*>(p);
	Eigen::Map<const EVector> nv(y, params->size);

	// Use these wrappers to pass the densities and overwrite the contents of dydt
	params->chemistry->evaluateFv(dydt, nv, *params->rateCoeffv, *params->kv);
	return GSL_SUCCESS; // maybe check for nan here
}

int ode_j(double /* unused t */, const double y[], double* dfdy, double dfdt[], void* p)
{
	auto* params = static_cast<struct ode_params*>(p);
	Eigen::Map<const EVector> nv(y, params->size);

	// The jacobian dfdy needs to be stored row-major for GSL.
	params->chemistry->evaluateJvv(dfdy, nv, *params->rateCoeffv, *params->Jkvv);

	// no explicit time dependence
	std::fill(dfdt, dfdt + params->size, 0);

	return GSL_SUCCESS;
}
} // namespace

void Chemistry::registerSpecies(const std::vector<std::string>& namev)
{
	_speciesIndex = SpeciesIndex(namev);
}

void Chemistry::addSpecies(const std::string& name) { _speciesIndex.addSpecies(name); }

void Chemistry::addReaction(const std::string& reactionName,
                            const std::vector<std::string>& reactantNamev,
                            const Array& reactantStoichv,
                            const std::vector<std::string>& productNamev,
                            const Array& productStoichv)
{
	Error::equalCheck("Lengths of list of species names and vector of coefficients",
	                  reactantNamev.size(), reactantStoichv.size());
	Error::equalCheck("Lengths of list of species names and vector of coefficients",
	                  productNamev.size(), productStoichv.size());

	// Give the reaction a number, and put its name in the map.
	_reactionIndexm.emplace(reactionName, _reactionIndexm.size());

	// Register the names and ratios of the species involved. Vectors with coefficients of
	// the right size will be created later.
	_reactionv.emplace_back(reactantNamev, reactantStoichv, productNamev, productStoichv);
}

void Chemistry::prepareCoefficients()
{
	_numSpecies = _speciesIndex.size();
	_rStoichvv = makeReactantStoichvv();
	_netStoichvv = makeProductStoichvv() - _rStoichvv;
	_numReactions = _rStoichvv.cols();
}

int Chemistry::reactionIndex(const std::string& reactionName) const
{
	return _reactionIndexm.at(reactionName);
}

EVector Chemistry::solveBalance(const EVector& rateCoeffv, const EVector& n0v,
                                double maxTime) const
{
	Error::equalCheck<int>("chemistry coefficients", _numReactions, _reactionv.size());
	Error::equalCheck<int>("chemistry coefficients", _numSpecies, _netStoichvv.rows());
	EVector result = solveTimeDep(rateCoeffv, n0v, maxTime);
	return result;
}

void Chemistry::evaluateFv(double* FvOutput, const EVector& nv, const EVector& rateCoeffv,
                           EVector& kv) const

{
	// This memory should be already allocated
	Eigen::Map<EVector> Fv(FvOutput, _numSpecies);

	// The total rate of each reaction, a.k.a. k(T) multiplied with the correct density
	// factor n_s ^ R_s,r (see notes).
	for (int r = 0; r < _numReactions; r++)
		kv(r) = reactionSpeed(nv, rateCoeffv, r);

	// Matrix multiplication between the net stoichiometry matrix (indexed on species,
	// reaction) and the rate vector (indexed on reaction).
	Fv = _netStoichvv * kv;
}

void Chemistry::evaluateJvv(double* JvvDataRowMajor, const EVector& nv,
                            const EVector& rateCoeffv, EMatrix& Jkvv) const

{
	// This memory should be already allocated. The elements are d f_i / d n_j, so the
	// jacobian should have (dimension of fv, dimension of nv = _numSpecies, _numSpecies).
	Eigen::Map<EMatrixRM> JFvv(JvvDataRowMajor, _numSpecies, _numSpecies);

	reactionSpeedJacobian(Jkvv, nv, rateCoeffv);

	// (d f_i / d_nj) = sum_r S_ir * (d k_r / d n_j)
	JFvv = _netStoichvv * Jkvv;
}

EMatrix Chemistry::makeReactantStoichvv() const
{
	EMatrix r(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
		r.col(j) = _speciesIndex.linearCombination(_reactionv[j]._rNamev,
		                                           _reactionv[j]._rCoeffv);
	return r;
}

EMatrix Chemistry::makeProductStoichvv() const
{
	EMatrix p(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
		p.col(j) = _speciesIndex.linearCombination(_reactionv[j]._pNamev,
		                                           _reactionv[j]._pCoeffv);
	return p;
}

EVector Chemistry::solveTimeDep(const EVector& rateCoeffv, const EVector& n0v,
                                double maxTime) const
{
	if (maxTime == 0)
		return n0v;

	bool toEquilibrium = maxTime < 0;

	// Parameters and workspace for ode_f. The latter will write to kTotalv sometimes, using
	// the pointer given to the struct below.
	EVector kv(_numReactions);
	EMatrix Jkvv(_numReactions, _numSpecies);
	struct ode_params params = {_numSpecies, &rateCoeffv, this, &kv, &Jkvv};

	gsl_odeiv2_system s{ode_f, ode_j, _numSpecies, &params};

	// 1 second initial step
	double ini_step = 1;
	double epsabs = 1;
	double epsrel = 1e-15;
	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&s, gsl_odeiv2_step_bsimp,
	                                                     ini_step, epsabs, epsrel);

	int max_steps = toEquilibrium ? 32 : 1;
	double t = 0;
	EVector nv = n0v;
	for (int i = 1; i <= max_steps; i++)
	{
		// If not going to equilibrium, we will take only one step, and the maxTime
		// argument will be used.
		double goalTime = maxTime;
		if (toEquilibrium)
		{
			// We need to advance time by enough to see a significant difference in the
			// slowest changing densities. If the chosen time scale is too short, then we
			// would wrongly conclude that some densities have reached equilibrium, while in
			// fact they just wouldn't have had enough time to change. Time scale =
			// max(density / rate of change). Don't forget abs because rate of change can be
			// negative of course.
			EArray fv(_numSpecies, 1);
			evaluateFv(fv.data(), nv, rateCoeffv, kv);
			double timeScale =
			                (fv != 0).select(nv.array() / fv.abs(), 0).maxCoeff();

			// Limit the time scale to roughly the age of the universe
			timeScale = std::min(timeScale, 5.5e17);

			// I found that using twice the value works slightly better, but there's a lot
			// of wiggle room of course
			goalTime = t + 2 * timeScale;
		}

		EVector previousNv = nv;
		// TODO: deal with possible "singular matrix" error from GSL. Usually appears
		// when integration time is too long, and no changes happen anymore (related to
		// precision?)
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

double Chemistry::reactionSpeed(const EVector& nv, const EVector& rateCoeffv, size_t r) const
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
	return densityProduct * rateCoeffv(r);
}

// JdensityProduct is indexed on (reaction r, d / d n_j)
void Chemistry::reactionSpeedJacobian(EMatrix& Jkvv, const EVector& nv,
                                      const EVector& rateCoeffv) const
{
	// Every column j is basically a reaction speed before it is multiplied with the rate
	// coefficient, derived with respect to one of the densities.
	Array densityPowers(_numSpecies);
	for (size_t r = 0; r < _numReactions; r++)
	{
		// precalculate these powers
		for (size_t s = 0; s < _numSpecies; s++)
		{
			double R = _rStoichvv(s, r);
			if (R)
				densityPowers[s] = pow(nv(s), R);
			else
				densityPowers[s] = 0;
		}

		// Now start calculating all the derivatives with respect to n_j. We cannot just
		// divide the product of the densities by n_j, because n_j can be zero.
		for (size_t j = 0; j < _numSpecies; j++)
		{
			double& Jrj = Jkvv(r, j);

			// if n_j not involved, just set to 0
			double Rj = _rStoichvv(j, r);
			if (Rj == 0)
			{
				Jrj = 0;
				continue;
			}
			// Otherwise, start with the rate coefficient for this reaction, and
			// then multiply with all the densities
			Jrj = rateCoeffv(r);

			// All species involved in the reaction except n_j
			for (size_t s = 0; s < _numSpecies; s++)
			{
				// if species is involved, multiply with n_s^Rs
				if (_rStoichvv(s, r) && s != j)
					Jrj *= densityPowers[s];
			}

			// derivative of the n_j^Rj factor: Rj n_j^(Rj - 1). Do not call pow if
			// the expornent is trivial (Rj == 1), or if Jrj is already zero.
			if (Rj != 1 && Jrj)
				Jrj = Rj * pow(nv(j), Rj - 1);
		}
	}
}
