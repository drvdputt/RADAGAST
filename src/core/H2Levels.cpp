#include "H2Levels.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GasStruct.hpp"
#include "H2FromFiles.hpp"
#include "LevelSolver.hpp"
#include "TemplatedUtils.hpp"

using namespace std;

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
