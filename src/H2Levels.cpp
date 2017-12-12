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

	// Initial guess (TODO: better initial guess with actual temperature? Maybe we don't even
	// need this).
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
		DEBUG(" lvl it " << counter << " norm " << previousNv.sum());
	}
	DEBUG("\nnv\n" << nv << std::endl);
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
	// used in chemical network (it will multiply by the density again.
	// TODO: need separate rates for H2g and H2*
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

	// Photon flux density in s-1 cm-2 Hz-1: 4pi I_nu / h nu
	Array photonFluxv =
	                Constant::FPI * specificIntensityv / Constant::PLANCK / frequencyv();

	// For each level
	for (size_t iLv = 0; iLv < numLv; iLv++)
	{
		// Integrate over the dissociation cross section TODO: optimize by only
		// integrating over the nonzero range. Best way to do this: provide raw access
		// to spectrum objects.
		Array sigmaTimesFlux{photonFluxv};
		for (size_t iNu = 0; iNu < frequencyv().size(); iNu++)
			sigmaTimesFlux *= _hff->directDissociationCrossSection(
			                frequencyv()[iNu], iLv);
		result(iLv) = TemplatedUtils::integrate<double, Array, Array>(frequencyv(), sigmaTimesFlux);
	}
	return result;
}

EVector H2Levels::spontaneousDissociationSinkv() const
{
	return EVector::Zero(_hff->numLv(), _hff->numLv());
}
