#include "HModel.hpp"

void HModel::solve(const GasStruct& gas, const Spectrum& specificIntensity) {}

Array HModel::emissivity(const Array eFrequencyv) const {}

Array HModel::opacityv(const Array oFrequencyv) const {}

double HModel::netHeating() const {}

Array HModel::twoPhotonEmissivityv(const Array& eFrequencyv) const
{
	std::array<size_t, 2> upper_lower = _hData->twoPhotonIndices();
	// This index can mean either the resolved level 2s, or the collapsed level 2
	size_t index2sOr2 = upper_lower[0];
	size_t index1s = upper_lower[1];

	const EVector& gv = _hData->gv();
	const EVector& ev = _hData->ev();

	// The population of the 2s level needs to be guessed when n=2 is collapsed. We can
	// check this by looking at the multiplicity of the level. Since 2s and 1s should both
	// have the the same multiplicity, we can do:
	bool collapsed = gv(index2sOr2) != gv(index1s);

	// If collapsed, assume the population of the 2s level to be 1/4 of the total n=2
	// population.
	double n2s = _levelSolution.nv()(index2sOr2);
	if (collapsed)
		n2s /= 4.;

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
