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
	EVector nv;
	if (gas._h2Levelv.size() == 0)
		nv = n * solveBoltzmanEquations(gas._T);
	else if (gas._h2Levelv.size() == numLv())
		nv = gas._h2Levelv;
	else
		Error::runtime("Wrong size for initial guess vector!");

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
                                  const Spectrum& specificIntensity) const
{
#ifdef STERNBERG2014
	// See 2014-Sternberg eq 3
	auto iv = specificIntensity.valuev();
	auto nuv = specificIntensity.frequencyv();

	// F0 = integral 912 to 1108 Angstrom of Fnu(= 4pi Inu) with Inu in cm-2 s-1 Hz sr-1
	Array photonFluxv = Constant::FPI * iv / nuv / Constant::PLANCK;
	constexpr double freqLWmin{Constant::LIGHT / 1108 / Constant::ANG_CM};
	constexpr double freqLWmax{Constant::LIGHT / 912 / Constant::ANG_CM};
	size_t iLWmin{TemplatedUtils::index(freqLWmin, nuv)};
	size_t iLWmax{TemplatedUtils::index(freqLWmax, nuv)};
	double F0 = TemplatedUtils::integrate<double>(nuv, photonFluxv, iLWmin, iLWmax);

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
		return dissociationSinkv(specificIntensity).dot(s.nv) / s.nv.sum();
	else
	{
		// We need to return something nonzero here, otherwise the chemistry will have
		// no dissociation coefficient, maxing out the H2.

		// Just pick the one for the ground state? Or maybe LTE? TODO: choose
		EVector lteRatios = solveBoltzmanEquations(s.T);
		return dissociationSinkv(specificIntensity).dot(lteRatios);
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

EVector H2Levels::dissociationSinkv(const Spectrum& specificIntensity) const
{
	return directDissociationSinkv(specificIntensity) + spontaneousDissociationSinkv();
}

EVector H2Levels::directDissociationSinkv(const Spectrum& specificIntensity) const
{
	size_t numLv{_hff->numLv()};
	EVector result{EVector::Zero(numLv)};

	// For each level that has cross section data
	for (size_t iLv : _levelsWithCrossSectionv)
	{
		// For each cross section
		for (const Spectrum& cs : _hff->directDissociationCrossSections(iLv))
		{
			// We will integrate over (part of, in case the input spectrum is not
			// wide enough) the grid for the cross section
			const Array& cs_nuv = cs.frequencyv();

			// Usable integration range
			double minFreq = max(cs.freqMin(), specificIntensity.freqMin());
			double maxFreq = min(cs.freqMax(), specificIntensity.freqMax());

			// Integration lower bound: Index right of the minimum frequency
			size_t iNuMin = TemplatedUtils::index(minFreq, cs_nuv);
			// Index left of the minimum frequency
			if (iNuMin > 0)
				iNuMin--;

			// Integration upper bound: Index right of the maximum frequency
			size_t iNuMax = TemplatedUtils::index(maxFreq, cs_nuv);

			// Integrand: flux * sigma
			// Start with cross section [cm-2]
			Array sigmaFv{cs.valuev()};
			for (size_t iNu = iNuMin; iNu <= iNuMax; iNu++)
			{
				// Multiply with photon flux density [s-1 cm-2 Hz-1]: F_nu = 4pi
				// I_nu / h nu. (As always constant factors are applied after
				// integrating.)
				double nu = cs_nuv[iNu];
				sigmaFv[iNu] *= specificIntensity.evaluate(nu) / nu;
			}
			// Integrate to total number of dissociations (s-1)
			result(iLv) += Constant::FPI / Constant::PLANCK *
			               TemplatedUtils::integrate<double>(cs_nuv, sigmaFv,
			                                                 iNuMin, iNuMax);
		}
	}
	return result;
}

EVector H2Levels::spontaneousDissociationSinkv() const
{
	return _hff->dissociationProbabilityv();
}
