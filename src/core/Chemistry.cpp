#include "Chemistry.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
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
        params->chemistry->evaluateFv(dydt, y, *params->rateCoeffv, *params->kv);
        return GSL_SUCCESS;  // maybe check for nan here
    }

    int ode_j(double /* unused t */, const double y[], double* dfdy, double dfdt[], void* p)
    {
        auto* params = static_cast<struct ode_params*>(p);
        params->chemistry->evaluateJvv(dfdy, y, *params->rateCoeffv, *params->Jkvv);

        // no explicit time dependence
        std::fill(dfdt, dfdt + params->size, 0);
        return GSL_SUCCESS;
    }
}  // namespace

void Chemistry::registerSpecies(const std::vector<std::string>& namev)
{
    _speciesIndex = SpeciesIndex(namev);
}

void Chemistry::addSpecies(const std::string& name)
{
    _speciesIndex.addSpecies(name);
}

void Chemistry::addReaction(const std::string& reactionName, const std::vector<std::string>& reactantNamev,
                            const std::vector<int>& reactantStoichv, const std::vector<std::string>& productNamev,
                            const std::vector<double>& productStoichv)
{
    Error::equalCheck("Lengths of list of species names and vector of coefficients", reactantNamev.size(),
                      reactantStoichv.size());
    Error::equalCheck("Lengths of list of species names and vector of coefficients", productNamev.size(),
                      productStoichv.size());

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
    _netStoichvv = makeProductStoichvv() - _rStoichvv.cast<double>();
    _numReactions = _rStoichvv.cols();
}

int Chemistry::reactionIndex(const std::string& reactionName) const
{
    return _reactionIndexm.at(reactionName);
}

EVector Chemistry::solveBalance(const EVector& rateCoeffv, const EVector& n0v, double maxTime) const
{
    Error::equalCheck<int>("chemistry coefficients", _numReactions, _reactionv.size());
    Error::equalCheck<int>("chemistry coefficients", _numSpecies, _netStoichvv.rows());
    EVector result = solveTimeDep(rateCoeffv, n0v, maxTime);
    return result;
}

void Chemistry::evaluateFv(double* FvOutput, const double* nv, const EVector& rateCoeffv, EVector& kv) const

{
    // The total rate of each reaction, a.k.a. k(T) multiplied with the correct density
    // factor n_s ^ R_s,r (see notes).
    for (int r = 0; r < _numReactions; r++) kv(r) = reactionSpeed(nv, rateCoeffv, r);

    // Matrix multiplication between the net stoichiometry matrix (indexed on species, reaction)
    // and the rate vector (indexed on reaction).

    // Using noalias() shaves roughly 30 percent off the total time spent on the chemistry. It
    // tells Eigen that it is safe to write to FvOutput while calculating the product, so no
    // temporary value is needed. This gets rid of an automatically generated resize() call.
    // Normally, Eigen decides the best option by itself, but here it can't because it doesn't know
    // at compile time what FvOuput points to, so it chooses the safest option by default.
    Eigen::Map<EVector>(FvOutput, _numSpecies).noalias() = _netStoichvv * kv;
}

void Chemistry::evaluateJvv(double* JvvDataRowMajor, const double* nv, const EVector& rateCoeffv, EMatrix& Jkvv) const

{
    // The derivative of the reaction speed vector d k_r / d n_j
    reactionSpeedJacobian(Jkvv, nv, rateCoeffv);

    // (d f_i / d_nj) = sum_r S_ir * (d k_r / d n_j). EMatrixRM has to be used because GSL stores
    // dfdy in row major order.
    Eigen::Map<EMatrixRM>(JvvDataRowMajor, _numSpecies, _numSpecies).noalias() = _netStoichvv * Jkvv;
}

EMatrix_int Chemistry::makeReactantStoichvv() const
{
    EMatrix_int R = EMatrix_int::Zero(_numSpecies, _reactionv.size());
    for (size_t r = 0; r < _reactionv.size(); r++)
    {
        for (size_t i = 0; i < _reactionv[r]._rNamev.size(); i++)
        {
            int s = _speciesIndex.index(_reactionv[r]._rNamev[i]);
            R(s, r) = _reactionv[r]._rCoeffv[i];
        }
    }
    return R;
}

EMatrix Chemistry::makeProductStoichvv() const
{
    EMatrix P = EMatrix::Zero(_numSpecies, _reactionv.size());
    for (size_t r = 0; r < _reactionv.size(); r++)
    {
        for (size_t i = 0; i < _reactionv[r]._pNamev.size(); i++)
        {
            int s = _speciesIndex.index(_reactionv[r]._pNamev[i]);
            P(s, r) = _reactionv[r]._pCoeffv[i];
        }
    }
    return P;
}

EVector Chemistry::solveTimeDep(const EVector& rateCoeffv, const EVector& n0v, double maxTime) const
{
    if (maxTime == 0) return n0v;

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
    gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&s, gsl_odeiv2_step_bsimp, ini_step, epsabs, epsrel);

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
            evaluateFv(fv.data(), nv.data(), rateCoeffv, kv);
            double timeScale = (fv != 0).select(nv.array() / fv.abs(), 0).maxCoeff();

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

double Chemistry::reactionSpeed(const double* nv, const EVector& rateCoeffv, size_t r) const
{
    double densityProduct = 1;
    for (size_t s = 0; s < _numSpecies; s++)
    {
        // If the species in involved in this reaction (stoich on left side > 0), calculate the
        // density to the power of its stoichiometry in the reaction.
        int Rs = _rStoichvv(s, r);
        if (Rs) densityProduct *= gsl_pow_int(nv[s], Rs);
    }
    return densityProduct * rateCoeffv(r);
}

void Chemistry::reactionSpeedJacobian(EMatrix& Jkvv, const double* nv, const EVector& rateCoeffv) const
{
    // Every column j is basically a reaction speed before it is multiplied with the rate
    // coefficient, derived with respect to one of the densities.
    Array densityPowers(_numSpecies);
    for (size_t r = 0; r < _numReactions; r++)
    {
        // precalculate these powers
        for (size_t s = 0; s < _numSpecies; s++)
        {
            int Rs = _rStoichvv(s, r);
            if (Rs)
                densityPowers[s] = gsl_pow_int(nv[s], Rs);
            else
                densityPowers[s] = 0;
        }

        // Now start calculating all the derivatives with respect to n_j. We cannot just divide the
        // product of the densities by n_j, because n_j can be zero.
        for (size_t j = 0; j < _numSpecies; j++)
        {
            double& Jrj = Jkvv(r, j);

            // if n_j not involved, just set to 0
            int Rj = _rStoichvv(j, r);
            if (Rj == 0)
            {
                Jrj = 0;
                continue;
            }
            // Otherwise, start with the rate coefficient for this reaction, and then multiply with
            // all the densities
            Jrj = rateCoeffv(r);

            // All species involved in the reaction except n_j
            for (size_t s = 0; s < _numSpecies; s++)
            {
                // if species is involved, multiply with n_s^Rs
                if (_rStoichvv(s, r) && s != j) Jrj *= densityPowers[s];
            }

            // derivative of the n_j^Rj factor: Rj n_j^(Rj - 1). Do not call pow if the exponent is
            // trivial (Rj == 1), or if Jrj is already zero.
            if (Rj != 1 && Jrj) Jrj *= Rj * gsl_pow_int(nv[j], Rj - 1);
        }
    }
}
