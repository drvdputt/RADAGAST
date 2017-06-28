#include "HydrogenLevels.h"
#include "Constants.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <iostream>
#include <vector>

HydrogenLevels::HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp)
                : NLevel(hdp.get()), _hdp(hdp)
{
}

Array HydrogenLevels::emissivityv(const Solution& s) const
{
	return lineEmissivityv(s) + twoPhotonEmissivityv(s);
}

Array HydrogenLevels::twoPhotonEmissivityv(const Solution& s) const
{
	std::array<int, 2> upper_lower = _hdp->twoPhotonIndices();

	// This index can mean either the resolved level 2s, or the collapsed level 2
	int index2sOr2 = upper_lower[0];
	int index1s = upper_lower[1];

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

	Array result(frequencyv().size());
	for (size_t iFreq = 0; frequencyv()[iFreq] < nu0; iFreq++)
	{
		double y = frequencyv()[iFreq] / nu0;
		double y1miny = y * (1 - y);
		double pow4y1miny_gam = pow(4 * y1miny, gam);
		double Py = C * (y1miny * (1 - pow4y1miny_gam) +
		                 alpha * pow(y1miny, beta) * pow4y1miny_gam);
		result[iFreq] = constFactor * y * Py;
	}
	return result;
}
