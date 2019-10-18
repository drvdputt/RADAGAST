#include "HydrogenLevels.hpp"
#include "Constants.hpp"
#include "GasStruct.hpp"
#include "HydrogenDataProvider.hpp"
#include "IonizationBalance.hpp"
#include "Options.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"

#include <vector>

using namespace std;

Array HydrogenLevels::emissivityv(const Solution& s, const Array& eFrequencyv) const
{
	return lineEmissivityv(s, eFrequencyv) + twoPhotonEmissivityv(s, eFrequencyv);
}

// namespace
// {
// double alpha(int n, int l, double T)
// {
// 	/* TODO: find better recombination coefficients. */

// 	/* for now use the hardcoded implementation, but this needs to change (is copy paste from
// 	   HydrogenHardcoded). */
// 	double T4 = T / 1.e4;
// 	double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4));
// 	double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
// 	double alpha2s = 2.34e-14 * pow(T4, -0.537 - 0.019 * log(T4));

// 	// 2015-Raga (A13)
// 	double t = log10(T4);
// 	vector<double> logAlpha3poly = {-13.3377, -0.7161, -0.1435, -0.0386, 0.0077};
// 	vector<double> logAlpha4poly = {-13.5225, -0.7928, -0.1749, -0.0412, 0.0154};
// 	vector<double> logAlpha5poly = {-13.6820, -0.8629, -0.1957, -0.0375, 0.0199};

// 	double alpha3 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha3poly));
// 	double alpha4 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha4poly));
// 	double alpha5 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha5poly));

// 	Array unresolvedAlphav({alphaGround, alpha2p + alpha2s, alpha3, alpha4, alpha5});

// 	if (n == 2 && l == 0)
// 		return alpha2s;
// 	else if (n == 2 && l == 1)
// 		return alpha2p;
// 	else if (n <= 5)
// 		return (2 * l + 1) * unresolvedAlphav[n - 1] / (2. * n * n);
// 	else
// 		return 0;
// }
// } // namespace
EVector HydrogenLevels::sourcev(const GasStruct& gas) const
{
	EVector result{EVector::Zero(_hdp->numLv())};

	// Now loop over all levels, and add the correct recombination coefficient. If a level
	// is collapsed, the alpha will be added to the same level multiple times.
	int nMax = _hdp->nMax();
	for (int n = 1; n <= nMax; n++)
	{
		for (int l = 0; l < n; l++)
		{
			size_t index = _hdp->indexOutput(n, l);
			result[index] += _rr->alpha(n, l, gas._T);
		}
	}

	if (Options::hlevels_topoff)
	{
		// The recombination coefficients should sum to the total one (the one used for the
		// ionization balance calculation). Therefore, we calculate here how much we still
		// need...
		double total_rr = Ionization::recombinationRateCoeff(gas._T);
		double topoff = total_rr - result.sum();
		if (topoff > 0)
		{
			// ...and add it to the top n-level
			for (int l = 0; l < nMax; l++)
			{
				size_t index = _hdp->indexOutput(nMax, l);
				result[index] += topoff / nMax / nMax * (2 * l + 1);
			}
		}
	}

	double ne = gas._speciesNv(SpeciesIndex::ine());
	double np = gas._speciesNv(SpeciesIndex::inp());
	return result * ne * np;
}

EVector HydrogenLevels::sinkv(const GasStruct& gas) const
{
	// TODO: ideally, this calculates the ionization rate from each level, using individual
	// ionization cross sections.

	/* The ionization rate calculation makes no distinction between the levels.  When
	   the upper level population is small, and its decay rate is large, the second term
	   doesn't really matter. Therefore, we choose the sink to be the same for each
	   level.  Moreover, total source = total sink so we want sink*n0 + sink*n1 = source
	   => sink = totalsource / n because n0/n + n1/n = 1. */
	double ne = gas._speciesNv(SpeciesIndex::ine());
	double np = gas._speciesNv(SpeciesIndex::inp());
	double nH = gas._speciesNv(SpeciesIndex::inH());
	double totalSource = ne * np * Ionization::recombinationRateCoeff(gas._T);
	double sink = totalSource / nH; // Sink rate per (atom per cm3)
	return EVector::Constant(numLv(), sink);
}
