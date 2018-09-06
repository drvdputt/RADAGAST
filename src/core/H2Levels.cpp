#include "H2Levels.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "GasStruct.h"
#include "H2FromFiles.h"
#include "TemplatedUtils.h"

using namespace std;

H2Levels::H2Levels(shared_ptr<const H2FromFiles> hff)
                : NLevel(hff, 2 * Constant::HMASS_CGS), _hff(hff)
{
	for (size_t i = 0; i < _hff->numLv(); i++)
		if (!_hff->directDissociationCrossSections(i).empty())
			_levelsWithCrossSectionv.emplace_back(i);
}

H2Levels::~H2Levels() = default;

Array H2Levels::opacityv(const Solution& s, const Array& oFrequencyv) const
{
	// Start with the line opacity
	Array totalOpv = lineOpacityv(s, oFrequencyv);

	// Then add the dissociation cross sections of each level
	for (size_t iLv : _levelsWithCrossSectionv)
	{
		const vector<Spectrum>& csv = _hff->directDissociationCrossSections(iLv);
		for (const Spectrum& cs : csv)
			totalOpv += cs.binned(oFrequencyv);
	}
	return totalOpv;
}

EVector H2Levels::solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
                                     const EVector& sourcev, const EVector& sinkv,
                                     int /* chooseConsvEq */, const GasStruct& gas) const
{
	// Initial guess
	EVector nv;
	if (gas._h2Levelv.size() == 0)
	{
		DEBUG("Using LTE as initial guess for H2" << endl);
		// nv = n * solveBoltzmanEquations(gas._T);
		// TODO: experiment with this instead of LTE as initial condition
		nv = EVector::Zero(numLv());
		nv(1) = n;
	}
	else if (gas._h2Levelv.size() == numLv())
		nv = gas._h2Levelv;
	else
		Error::runtime("Wrong size for initial guess vector!");

	// This should stay constant during the calculation
	const EMatrix Mvv = netTransitionRate(BPvv, Cvv);

	// Fractional destruction rate (in s-1) stays constant when populations are adjusted
	// (sum over the row indices by doing a colwise reduction. Need to transpose because
	// summing each column gives a row vector, while sinkv is column one.
	EVector fracDestructionRatev = Mvv.colwise().sum().transpose() + sinkv;

	// Get the indices that cover the electronic ground state
	auto indicesX = _hff->indicesX();
	auto indicesExcited = _hff->indicesExcited();

	// The algorithm apparently works better when the iterations happen from high to low.
	// While it is not guaranteed that the indices given by the above vector are actually
	// listed from low to high energies, it is so in practice (because of the way the files
	// that are read in are organized. Therefore, we reverse the vectors here, and hope for
	// the best, as currently I'm experiencing oscillations under strong radiation fields.
	reverse(begin(indicesX), end(indicesX));
	reverse(begin(indicesExcited), end(indicesExcited));

	// Iterate until converged
	double maxDeltaX = 1.e-4 * n;
	double maxDeltaAll = 1.e-5 * n;
	double maxFracDelta = 1.e-2;
	size_t counter{0};
	double xFrac = 0, eFrac = 0;
	while (true)
	{
		counter++;

		// The previous nv, for convergence checking
		EVector previousNv = nv;
		EVector deltav = EVector::Zero(nv.size());

		// Sweep over ground state
		xFrac = 0;
		for (auto i : indicesX)
		{
			// Sum Mij nj, with j running over all other levels. TODO: consider
			// making the assumption that the X states are contigous in the eigen
			// index space. In that case, I can simply do a matrix operation here,
			// over a slice, leading to possible optimization.
			double creationRate = sourcev(i) + Mvv.row(i) * nv;
			nv(i) = creationRate <= 0 ? 0 : creationRate / fracDestructionRatev(i);
			xFrac += nv(i);
		}
		/* Renormalize because the algorithm has no sum rule, */
		nv *= n / nv.sum();
		xFrac /= n;

		// We will cut this iteration short as long as the ground state has not
		// converged.
		bool loopBack = false;
		for (auto i : indicesX)
		{
			// TODO: vectorize this (maybe work with a contiguous 'range' of X
			// indices, instead of a list
			deltav(i) = abs(nv(i) - previousNv(i));
			if (!loopBack && deltav(i) > maxDeltaX)
			{
				// DEBUG("index " << i << " delta " << deltav(i) / maxDeltaX << endl);
				// Do not break here. We might need the whole deltav further down
				loopBack = true;
			}
		}
		if (loopBack)
		{
			if (!(counter % 1000))
			{
				DEBUG("Solving h2... " << counter << "iterations" << endl);
				DEBUG("X fraction = " << xFrac << endl);
				// DEBUG("h2Levelv = \n" << nv << endl);
			}
			continue;
		}

		// If the ground state has more or less converged, we will also start sweeping
		// over the other states
		eFrac = 0;
		for (auto i : indicesExcited)
		{
			double creationRate = sourcev(i);
			// Only sum over the ground state here
			for (auto j : indicesX)
				creationRate += Mvv(i, j) * nv(j);
			nv(i) = creationRate <= 0 ? 0 : creationRate / fracDestructionRatev(i);
			eFrac += nv(i);
		}
		nv *= n / nv.sum();
		eFrac /= n;

		// Overall convergence check

		// Absolute change
		for (auto i : indicesExcited)
			// Again, this would best be vectorized
			deltav(i) = nv(i) - previousNv(i);

		bool thresconv = (deltav.array() < maxDeltaAll).all();

		// Relative change
		bool fracconv = true;
		for (int i = 0; i < nv.size() && fracconv; i++)
		{
			double df = 0;
			// finite/0 means not converged
			if (previousNv(i) <= 0 && nv(i) > 0)
				df = 2 * maxFracDelta;
			else
				df = abs(nv(i) / previousNv(i) - 1.);

			if (df > maxFracDelta)
				fracconv = false;
		}

		if (!(counter % 1000))
		{
			DEBUG("Solving h2... " << counter << "iterations" << endl);
			DEBUG("thresconv = " << thresconv << " fracconv = " << fracconv
			                     << endl);
		}

		bool converged = thresconv && fracconv;
		const int max_iterations = 10000;
		if (converged || counter > max_iterations)
			break;
	}
	if (nv.array().isNaN().any())
		Error::runtime("nan in H2 level solution");
	DEBUG("Solved H2 in " << counter << " iterations" << endl);
	DEBUG("Excited fraction = " << eFrac << endl);
	// DEBUG("h2Levelv = \n" << nv << endl);
	// EVector explicitNv = NLevel::solveRateEquations(n, BPvv, Cvv, sourcev, sinkv,
	// chooseConsvEq, gas);
	return nv;
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

double H2Levels::dissociationCooling(const Solution&) const { return 0.0; }

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
