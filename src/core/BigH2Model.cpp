#include "BigH2Model.hpp"
#include "GasStruct.hpp"
#include "LevelSolver.hpp"

void BigH2Model::solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
                       double h2form)
{
	_n = n;
	_levelSolution.setT(gas._T);

	EMatrix Cvv;
	EMatrix Tvv = _h2Data->totalTransitionRatesvv(specificIntensity, gas, &Cvv);
	_levelSolution.setCvv(Cvv);

	// TODO: Use better formation pumping recipe
	EVector sourcev = EVector::Zero(_h2Data->numLv());
	sourcev.head(_h2Data->startOfExcitedIndices()) = _h2Data->formationDistribution();
	sourcev *= h2form / sourcev.sum();

	EVector sinkv = dissociationSinkv(specificIntensity);

	if (_levelSolution.nv().size() == 0)
	{
		EVector initialGuessv = EVector::Zero(_h2Data->numLv());

		// Use LTE for the X levels, and 0 for the rest
		int endX = _h2Data->startOfExcitedIndices();
		initialGuessv.head(endX) = LevelSolver::statisticalEquilibrium_boltzman(
		                n, gas._T, _h2Data->ev().head(endX), _h2Data->gv().head(endX));

		DEBUG("Using LTE as initial guess for H2" << std::endl);
		_levelSolution.setNv(initialGuessv);
	}

	int fullyConnectedCutoff = _h2Data->startOfExcitedIndices();
	EVector newNv = LevelSolver::statisticalEquilibrium_iterative(
	                n, Tvv, sourcev, sinkv, _levelSolution.nv(), fullyConnectedCutoff);
	_levelSolution.setNv(newNv);
}

double BigH2Model::dissociationRate(const Spectrum& specificIntensity) const
{
#ifdef STERNBERG2014
// TODO: move this to simple h2 model
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

	EVector popFracv;
	if (_n > 0)
		popFracv = _levelSolution.nv() / _n;
	else
		// We need to return something nonzero here, otherwise the chemistry will have
		// no dissociation coefficient, which can be troublesomec.
		popFracv = _h2Data->solveBoltzmanEquations(_levelSolution.t());

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

double BigH2Model::orthoPara() const
{
	// TODO
	return .25;
}

Array BigH2Model::emissivityv(const Array& eFrequencyv) const
{
	return _levelSolution.emissivity(eFrequencyv);
}

Array BigH2Model::opacityv(const Array& oFrequencyv) const
{
	// Start with the line opacity
	Array totalOpv = _levelSolution.opacityv(oFrequencyv);

	// Then add the dissociation cross sections of each level
	for (size_t iLv : _h2Data->levelsWithCrossSectionv())
	{
		const std::vector<Spectrum>& csv =
		                _h2Data->directDissociationCrossSections(iLv);
		for (const Spectrum& cs : csv)
			totalOpv += cs.binned(oFrequencyv);
	}
	return totalOpv;
}
