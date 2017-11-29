#include "H2Levels.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "H2FromFiles.h"
#include "TemplatedUtils.h"

H2Levels::H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv)
                : NLevel(hff, frequencyv), _hff(hff)
{
}

H2Levels::~H2Levels() = default;

/** Mvv is Mij in my notes. */
double evaluateSinglePopulation(size_t i, const EVector& nv, const EMatrix& Mvv,
                                const EVector& sourcev, const EVector& sinkv)
{
	// TODO: different implementation for X vs excited levels
	// creationRate - ni * destructionFraction = 0
	// ==> ni = creationRate / destructionFraction

	// FIXME: assume that we only have X states (this is the case at the moment (28/11/2017))) for excited states the summation boundaries are different

	// sum_j M_ij n_j + f_i
	double creationRate = (Mvv.row(i) * nv).sum() + sourcev(i);
	if (creationRate <= 0)
		return 0;

	// sum_j M_ji + d_i
	double destructionFraction = Mvv.col(i).sum() + sinkv(i);

	return creationRate / destructionFraction;
}

EVector H2Levels::solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
                                     const EVector& sourcev, const EVector& sinkv,
                                     int chooseConsvEq) const
{
#define USE_ITERATION_METHOD
#ifdef USE_ITERATION_METHOD
	// This should stay constant during the calculation
	const EMatrix Mvv = netTransitionRate(BPvv, Cvv);

	// Initial guess (TODO: better initial guess with actual temperature? Maybe we don't even
	// need this).
	EVector nv = n * solveBoltzmanEquations(100);

	// Iterate until converged
	bool converged = false;
	size_t counter{0};
	while (!converged)
	{
		EVector previousNv = nv;
		for (size_t i = 0; i < _hff->numLv(); i++)
			nv(i) = evaluateSinglePopulation(i, nv, Mvv, sourcev, sinkv);

		/* Renormalize; the algorithm has no sum rule, */
		nv *= n / nv.sum();

		EVector deltaNv = (nv - previousNv).array().abs();
		converged = deltaNv.maxCoeff() < 1e-6 * n;
		counter++;
		DEBUG(" lvl it " << counter << " norm " << previousNv.sum());
	}
	DEBUG("\nnv\n" << nv << std::endl);
	return nv;
#else
	return NLevel::solveRateEquations(n, BPvv, Cvv, sourcev, sinkv, chooseConsvEq);
#endif
}

double H2Levels::dissociationRate(const NLevel::Solution& s, const Array& specificIntensityv) const
{
	// See 2014-Sternberg eq 3
	// F0 = integral 912 to 1108 Angstrom of Fnu(= 4pi Inu) with Inu in cm-2 s-1 Hz sr-1
	Array photonFluxv = Constant::FPI * specificIntensityv / frequencyv() / Constant::PLANCK;
	constexpr double freqLWmin{Constant::LIGHT / 1108 / Constant::ANG_CM};
	constexpr double freqLWmax{Constant::LIGHT / 912 / Constant::ANG_CM};
	size_t iLWmin{TemplatedUtils::index(freqLWmin, frequencyv())};
	size_t iLWmax{TemplatedUtils::index(freqLWmax, frequencyv())};
	double F0 = TemplatedUtils::integrate<double>(frequencyv(), photonFluxv, iLWmin, iLWmax);

	// eq 4 and 5
	double Iuv{F0 / 2.07e7};
	double result{5.8e-11 * Iuv};

	return result;
}

double H2Levels::dissociationHeating(const NLevel::Solution& s) const
{
	// TODO
	return 0;
}
