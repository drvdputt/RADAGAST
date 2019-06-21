#include "H2Levels.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "GasStruct.h"
#include "H2FromFiles.h"
#include "LevelSolver.h"
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

NLevel::Solution H2Levels::customSolution(double n, const GasStruct& gas,
                                          const Spectrum& specificIntensity,
                                          double h2form) const
{
	NLevel::Solution s;
	s.n = n;
	s.T = gas._T;
	EMatrix Tvv = totalTransitionRatesvv(specificIntensity, gas, &s.cvv);

	// TODO: Use better formation pumping recipe
	EVector sourcev = EVector::Zero(numLv());
	sourcev.head(_hff->startOfExcitedIndices()) = _hff->formationDistribution();
	sourcev *= h2form / sourcev.sum();

	EVector sinkv = dissociationSinkv(specificIntensity);

	EVector initialGuessv;
	if (gas._h2Levelv.size() == 0)
	{
		initialGuessv = EVector::Zero(numLv());

		// Use LTE for the X levels, and 0 for the rest
		int endX = _hff->startOfExcitedIndices();
		initialGuessv.head(endX) = LevelSolver::statisticalEquilibrium_boltzman(
		                n, gas._T, ev().head(endX), gv().head(endX));

		DEBUG("Using LTE as initial guess for H2" << endl);
	}
	else if (gas._h2Levelv.size() == numLv())
		initialGuessv = gas._h2Levelv;
	else
		Error::runtime("The _h2levelv in the gas struct should be empty, or be of the "
		               "right size already.");

	int fullyConnectedCutoff = _hff->startOfExcitedIndices();
	s.nv = LevelSolver::statisticalEquilibrium_iterative(
	                n, Tvv, sourcev, sinkv, initialGuessv, fullyConnectedCutoff);
	// s.nv = LevelSolver::statisticalEquilibrium_iterative(n, Tvv, sourcev, sinkv,
	//                                                      initialGuessv);
	return s;
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
	EVector directv = directDissociationIntegralv(specificIntensity);
	EVector solomonv = spontaneousDissociationSinkv();

	double nH2 = s.nv.sum();

	EVector popFracv;
	if (nH2 > 0)
		popFracv = s.nv / nH2;
	else
		// We need to return something nonzero here, otherwise the chemistry will have
		// no dissociation coefficient, which can be troublesomec.
		popFracv = solveBoltzmanEquations(s.T);

	// Dot product = total rate [cm-3 s-1]. Divide by total to get [s-1] rate, which
	// can be used in chemical network (it will multiply by the density again. TODO:
	// need separate rates for H2g and H2*
	double directFractional = directv.dot(popFracv);
	double solomonFractional = solomonv.dot(popFracv);
	DEBUG("Dissociation: direct rate:" << directFractional
	                                   << " solomon rate: " << solomonFractional << '\n');
	return directFractional + solomonFractional;
#endif
}

double H2Levels::dissociationHeating(const Solution& s, const Spectrum& specificIntensity) const
{
	// Fraction that dissociates per second * kinetic energy per dissociation * density of
	// level population = heating power density
	EArray p = _hff->dissociationProbabilityv().array();
	EArray k = _hff->dissociationKineticEnergyv().array();
	double solomonHeat = s.nv.dot(EVector{p * k});

	EVector directHeatv = directDissociationIntegralv(specificIntensity, true);
	double directHeat = s.nv.dot(directHeatv);

	DEBUG("Dissociation heat: direct heat: " << directHeat
	                                         << " solomon heat:" << solomonHeat << '\n');
	return solomonHeat + directHeat;
}

double H2Levels::dissociationCooling(const Solution&) const { return 0.0; }

EVector H2Levels::dissociationSinkv(const Spectrum& specificIntensity) const
{
	EVector directv = directDissociationIntegralv(specificIntensity);
	EVector solomonv = spontaneousDissociationSinkv();
	return directv + solomonv;
}

EVector H2Levels::directDissociationIntegralv(const Spectrum& specificIntensity,
                                              bool heatRate) const
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

				// If we are calculating the heating rate (erg s-1) instead of
				// the number rate (s-1), multiply with the energy minus the
				// threshold (i.e. the lowest frequency of the grid for that
				// specific cross section)
				if (heatRate)
					sigmaFv[iNu] *= Constant::PLANCK * (nu - cs_nuv[0]);
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
