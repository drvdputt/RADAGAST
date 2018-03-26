#include "H2Levels.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "GasStruct.h"
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
                                     int chooseConsvEq, const GasStruct& gas) const
{
#define USE_ITERATION_METHOD
#ifdef USE_ITERATION_METHOD
	// Initial guess
	EVector nv = n * solveBoltzmanEquations(gas._T);

	// This should stay constant during the calculation
	const EMatrix Mvv = netTransitionRate(BPvv, Cvv);

	// Fractional destruction rate (in s-1) stays constant when populations are adjusted
	// (sum over the row indices by doing a colwise reduction. Need to transpose because
	// summing each column gives a row vector, while sinkv is column one.
	EArray fracDestructionRatev = Mvv.colwise().sum().transpose() + sinkv;

	// Get the indices that cover the electronic ground state
	const auto& indicesX = _hff->indicesX();
	const auto& indicesExcited = _hff->indicesExcited();

	// Iterate until converged
	bool converged = false;
	size_t counter{0};
	while (!converged)
	{
		// The previous nv, for convergence checking
		EVector previousNv = nv;

		// Sweep over ground state
		for (auto i : indicesX)
		{
			// Sum Mij nj, with j running over all other levels. TODO: consider
			// making the assumption that the X states are contigous in the eigen
			// index space. In that case, I can simply do a matrix operation here,
			// over a slice, leading to possible optimization.
			double creationRate = Mvv.row(i) * nv + sourcev(i);
			if (creationRate <= 0)
				nv(i) = 0;
			else
				nv(i) = creationRate / fracDestructionRatev(i);
		}

		// Sweep over the other states
		for (auto i : indicesExcited)
		{
			double creationRate = sourcev(i);
			// Only sum over the ground state here
			for (auto j : indicesX)
				creationRate += Mvv(i, j) * nv(j);
			if (creationRate <= 0)
				nv(i) = 0;
			else
				nv(i) = creationRate / fracDestructionRatev(i);
		}
		// TODO: if this is still too slow, try the following optimizations:

		// 1. First put all the creation rates in a vector, then try the division (this
		// makes dealing with zeros a little harder though)

		/* Renormalize because the algorithm has no sum rule, */
		nv *= n / nv.sum();

		EVector deltaNv = (nv - previousNv).array().abs();
		converged = deltaNv.maxCoeff() < 1e-6 * n;
		counter++;
		if (!(counter % 1000))
			DEBUG("Solving h2... " << counter << "iterations" << endl);
	}
	if (nv.array().isNaN().any())
		Error::runtime("nan in H2 level solution");
	DEBUG("Solved H2 in " << counter << " iterations" << endl);
	// DEBUG("\nnv\n" << nv << endl);
	return nv;
#else
	return NLevel::solveRateEquations(n, BPvv, Cvv, sourcev, sinkv, chooseConsvEq, gas);
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
	double nH2 = s.nv.sum();
	if (nH2 > 0)
		return dissociationSinkv(specificIntensityv).dot(s.nv) / s.nv.sum();
	else
	{
		// We need to return something nonzero here, otherwise the chemistry will have
		// no dissociation coefficient, maxing out the H2.

		// Just pick the one for the ground state? Or maybe LTE? TODO: choose
		EVector lteRatios = solveBoltzmanEquations(s.T);
		return dissociationSinkv(specificIntensityv).dot(lteRatios);
	}
#endif
}

double H2Levels::dissociationHeating(const Solution& s) const
{
	// Fraction that dissociates per second * kinetic energy per dissociation * density of
	// level population = heating power density
	EArray p = _hff->dissociationProbabilityv().array();
	EArray k = _hff->dissociationKineticEnergyv().array();
	return s.nv.dot(EVector{p * k});
}

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
				// Convert to photon flux density in s-1 cm-2 Hz-1: F_nu = 4pi
				// I_nu / h nu
				sigmaFv[j] = Constant::FPI * sigmaFv[j] / Constant::PLANCK /
				             nuv[j];
				// Convert to dissociation count (s-1 Hz-1)
				sigmaFv[j] *= cs.evaluate(nuv[j]);
			}
			// Integrate to total number of dissociations (s-1)
			result(iLv) += TemplatedUtils::integrate<double>(nuv, sigmaFv);
		}
	}
	return result;
}

EVector H2Levels::spontaneousDissociationSinkv() const
{
	return _hff->dissociationProbabilityv();
}
