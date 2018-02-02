#include "H2Levels.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "H2FromFiles.h"
#include "TemplatedUtils.h"

using namespace std;

H2Levels::H2Levels(shared_ptr<const H2FromFiles> hff, const Array& frequencyv)
                : NLevel(hff, frequencyv, 2 * Constant::HMASS_CGS), _hff(hff)
{
	for (size_t i = 0; i < _hff->numLv(); i++)
		if (!_hff->directDissociationCrossSections(i).empty())
			_levelsWithCrossSectionv.emplace_back(i);
}

H2Levels::~H2Levels() = default;

Array H2Levels::opacityv(const Solution& s) const
{
	// Start with the line opacity
	Array totalOpv = lineOpacityv(s);

	// Then add the dissociation cross section of each level, for each frequency
	const Array& freqv = frequencyv();
	for (size_t iLv = 0; iLv < _hff->numLv(); iLv++)
		for (size_t iNu = 0; iNu < freqv.size(); iNu++)
			totalOpv[iNu] += s.nv(iLv) *
			                 _hff->directDissociationCrossSection(freqv[iNu], iLv);
	return totalOpv;
}

EVector H2Levels::solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
                                     const EVector& sourcev, const EVector& sinkv,
                                     int chooseConsvEq) const
{
#define USE_ITERATION_METHOD
#ifdef USE_ITERATION_METHOD
	size_t numLv = _hff->numLv();

	// This should stay constant during the calculation
	const EMatrix Mvv = netTransitionRate(BPvv, Cvv);

	// Initial guess (TODO: better initial guess with actual temperature? Maybe we don't
	// even need this).
	EVector nv = n * solveBoltzmanEquations(100);

	// Destruction rate (in s-1) stays constant when populations are adjusted
	Array destructionRatev(numLv);
	for (size_t i = 0; i < numLv; i++)
		destructionRatev[i] = Mvv.col(i).sum() + sinkv(i);

	// Iterate until converged
	bool converged = false;
	size_t counter{0};
	while (!converged)
	{
		// The previous nv, for convergence checking
		EVector previousNv = nv;

		/* Do a 'sweep' over all the populations. It is important that this happens one
		   by one, and not as a single vector operation (hence the word sweep). */
		for (size_t i = 0; i < numLv; i++)
		{
			double creationRate = (Mvv.row(i) * nv).sum() + sourcev(i);
			if (creationRate <= 0)
				nv(i) = 0;
			else
				nv(i) = creationRate / destructionRatev[i];
		}

		/* Renormalize because the algorithm has no sum rule, */
		nv *= n / nv.sum();

		EVector deltaNv = (nv - previousNv).array().abs();
		converged = deltaNv.maxCoeff() < 1e-6 * n;
		counter++;
	}
	DEBUG("Solved H2 in " << counter << " iterations" << endl);
	// DEBUG("\nnv\n" << nv << endl);
	return nv;
#else
	return NLevel::solveRateEquations(n, BPvv, Cvv, sourcev, sinkv, chooseConsvEq);
#endif
}

double H2Levels::dissociationRate(const NLevel::Solution& s,
                                  const Array& specificIntensityv) const
{
#ifdef STERNBERG2014
	// See 2014-Sternberg eq 3

	// F0 = integral 912 to 1108 Angstrom of Fnu(= 4pi Inu) with Inu in cm-2 s-1 Hz sr-1
	Array photonFluxv =
	                Constant::FPI * specificIntensityv / frequencyv() / Constant::PLANCK;
	constexpr double freqLWmin{Constant::LIGHT / 1108 / Constant::ANG_CM};
	constexpr double freqLWmax{Constant::LIGHT / 912 / Constant::ANG_CM};
	size_t iLWmin{TemplatedUtils::index(freqLWmin, frequencyv())};
	size_t iLWmax{TemplatedUtils::index(freqLWmax, frequencyv())};
	double F0 = TemplatedUtils::integrate<double>(frequencyv(), photonFluxv, iLWmin,
	                                              iLWmax);

	// eq 4 and 5
	double Iuv{F0 / 2.07e7};
	double result{5.8e-11 * Iuv};
	return result;
#else
	// Dot product = total rate [cm-3 s-1]. Divide by total to get [s-1] rate, which can be
	// used in chemical network (it will multiply by the density again. TODO: need separate
	// rates for H2g and H2*
	return dissociationSinkv(specificIntensityv).dot(s.nv) / s.nv.sum();
#endif
}

double H2Levels::dissociationHeating(const Solution& s) const { return 0.0; }

double H2Levels::dissociationCooling(const Solution& s) const { return 0.0; }

EVector H2Levels::dissociationSinkv(const Array& specificIntensityv) const
{
	return directDissociationSinkv(specificIntensityv) + spontaneousDissociationSinkv();
}

EVector H2Levels::directDissociationSinkv(const Array& specificIntensityv) const
{
	size_t numLv{_hff->numLv()};
	EVector result{EVector::Zero(numLv)};

	// For each level that has cross section data
	for (size_t iLv : _levelsWithCrossSectionv)
	{
		// For each cross section
		for (const Spectrum& cs : _hff->directDissociationCrossSections(iLv))
		{
			// Integration lower bound: Index right of the minimum frequency
			size_t iNuMin = TemplatedUtils::index(cs.freqMin(), frequencyv());
			// Index left of the minimum frequency
			if (iNuMin > 0)
				iNuMin--;

			// Integration upper bound: Index right of the maximum frequency
			size_t iNuMax = TemplatedUtils::index(cs.freqMax(), frequencyv());

			// Integration points
			vector<double> nuv(begin(frequencyv()) + iNuMin,
			                   begin(frequencyv()) + iNuMax);

			// Integrand: flux * sigma
			// Start with specific intensity
			vector<double> sigmaFv(begin(specificIntensityv) + iNuMin,
			                       begin(specificIntensityv) + iNuMax);
			for (size_t j = 0; j < nuv.size(); j++)
			{
				// Convert to photon flux density in s-1 cm-2 Hz-1: F_nu = 4pi I_nu / h nu
				sigmaFv[j] = Constant::FPI * sigmaFv[j] / Constant::PLANCK / nuv[j];
				// Convert to dissociation count (s-1 Hz-1)
				sigmaFv[j] *= cs.evaluate(nuv[j]);
			}
			// Integrate to total number of dissociations (s-1)
			result(iLv) += TemplatedUtils::integrate<double>(nuv, sigmaFv);
		}
	}
	return result;
}

EVector H2Levels::spontaneousDissociationSinkv() const { return EVector::Zero(_hff->numLv()); }
