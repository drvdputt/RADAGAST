#include "HydrogenLevels.h"
#include "Constants.h"
#include "GasStruct.h"
#include "HydrogenDataProvider.h"
#include "IonizationBalance.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"

#include <vector>

using namespace std;

HydrogenLevels::HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp)
                : NLevel(hdp, Constant::HMASS_CGS), _hdp(hdp), _rr(make_unique<HydrogenADF48>())
{
}

HydrogenLevels::~HydrogenLevels() = default;

Array HydrogenLevels::emissivityv(const Solution& s, const Array& eFrequencyv) const
{
	return lineEmissivityv(s, eFrequencyv) + twoPhotonEmissivityv(s, eFrequencyv);
}

namespace
{
double alpha(int n, int l, double T)
{
	/* TODO: find better recombination coefficients. */

	/* for now use the hardcoded implementation, but this needs to change (is copy paste from
	   HydrogenHardcoded). */
	double T4 = T / 1.e4;
	double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4));
	double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
	double alpha2s = 2.34e-14 * pow(T4, -0.537 - 0.019 * log(T4));

	// 2015-Raga (A13)
	double t = log10(T4);
	vector<double> logAlpha3poly = {-13.3377, -0.7161, -0.1435, -0.0386, 0.0077};
	vector<double> logAlpha4poly = {-13.5225, -0.7928, -0.1749, -0.0412, 0.0154};
	vector<double> logAlpha5poly = {-13.6820, -0.8629, -0.1957, -0.0375, 0.0199};

	double alpha3 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha3poly));
	double alpha4 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha4poly));
	double alpha5 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha5poly));

	Array unresolvedAlphav({alphaGround, alpha2p + alpha2s, alpha3, alpha4, alpha5});

	if (n == 2 && l == 0)
		return alpha2s;
	else if (n == 2 && l == 1)
		return alpha2p;
	else if (n <= 5)
		return (2 * l + 1) * unresolvedAlphav[n - 1] / (2. * n * n);
	else
		return 0;
}
} // namespace

EVector HydrogenLevels::sourcev(const GasStruct& gas) const
{
	EVector result{EVector::Zero(_hdp->numLv())};

	// Now loop over all levels, and add the correct recombination coefficient. If a level
	// is collapsed, the alpha will be added to the same level multiple times.
	for (int n = 1; n <= _hdp->nMax(); n++)
	{
		for (int l = 0; l < n; l++)
		{
			size_t index = _hdp->indexOutput(n, l);
			result[index] += _rr->alpha(n, l, gas._T);
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

Array HydrogenLevels::twoPhotonEmissivityv(const Solution& s, const Array& eFrequencyv) const
{
	std::array<size_t, 2> upper_lower = _hdp->twoPhotonIndices();

	// This index can mean either the resolved level 2s, or the collapsed level 2
	size_t index2sOr2 = upper_lower[0];
	size_t index1s = upper_lower[1];

	/* The population of the 2s level needs to be guessed when n=2. We can check this by looking
	   at the multiplicity of the level. Since 2s and 1s should both have the the same
	   multiplicity, we can do: */
	bool collapsed = gv(index2sOr2) != gv(index1s);
	/* If collapsed, assume the population of the 2s level to be 1/4 of the total n=2
	   population. */
	double n2s = collapsed ? s.nv(index2sOr2) / 4. : s.nv(index2sOr2);

	// 1984-Nussbaumer
	// constant factor in eq 3
	double constFactor = Constant::PLANCK / Constant::FPI * n2s;
	double nu0 = (ev(index2sOr2) - ev(index1s)) / Constant::PLANCK;

	// Parameters for eq 2
	const double C = 202.0; // s-1
	const double alpha = .88;
	const double beta = 1.53;
	const double gam = .8;

	Array result(eFrequencyv.size());
	for (size_t iFreq = 0; eFrequencyv[iFreq] < nu0 && iFreq < eFrequencyv.size(); iFreq++)
	{
		double y = eFrequencyv[iFreq] / nu0;
		double y1miny = y * (1 - y);
		double pow4y1miny_gam = pow(4 * y1miny, gam);
		double Py = C * (y1miny * (1 - pow4y1miny_gam) +
		                 alpha * pow(y1miny, beta) * pow4y1miny_gam);
		result[iFreq] = constFactor * y * Py;
	}
	return result;
}
