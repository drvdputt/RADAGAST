#ifndef CORE_HMODEL_HPP
#define CORE_HMODEL_HPP

#include "HData.hpp"
#include "LevelSolution.hpp"

class HModel
{
public:
	/** Pass a pointer to the HData object at construction, so that the constant H data and
	    functions can be accessed. */
	HModel(const HData* hData) : _hData{hData}, _levelSolution(_hData) {}

	/** Solve the H levels, and store them in this object. */
	void solve(const GasStruct& gas, const Spectrum& specificIntensity);

	/** This function returns the line emission spectrum + the continuum emitted by the
	    2s-1s two-photon process. */
	Array emissivityv(const Array eFrequencyv) const;

	/** This function returns the line opacity. TODO: maybe include bound-free cross section
	    that depends on the level populations. */
	Array opacityv(const Array oFrequencyv) const;

	/** From the level populations, calculate the net heating-cooling balance by
	    (de-)excitation */
	double netHeating() const;

	/** Read acces to the level solution, should one wish to inspect the solution in more
	    detail */
	const LevelSolution* levelSolution() const { return &_levelSolution; }

private:
	/** Returns a vector containing the source terms for the equilibrium equations, such as the
	    partial recombination rates into each level. Note that this function is separated from
	    the ionization balance calculation, as there only the total recombination rate matters.
	    In this case, this function returns the partial recombination rate into each level,
	    interpolated from some formulae found in 2015-Raga (for levels 3-5) , and Draine's book
	    (for levels 1, 2). The l-resolved recombination rates are just weighed by the degeneracy
	    of the level, for levels 3, 4 and 5. TODO: Use better data here. [cm-3 s-1] */
	EVector sourcev(const GasStruct& gas) const;

	/** Produces the sink term to be used by the equilibrium equations. In this case, hydrogen
	    disappears from the level populations because it's being ionized. In the current
	    implementation, all the ionization is assumed to be drawn equally from all levels. TODO:
	    add the effects of H2 formation in here?. Take care of this using actual ionization
	    cross sections? [s-1] */
	EVector sinkv(const GasStruct& gas) const;

	/** This function calculates the two-photon continuum using Nussbaumer \& Smutz (1984). */
	Array twoPhotonEmissivityv(const Array& eFrequencyv) const;

	const HData* _hData;
	LevelSolution _levelSolution;
};

#endif // CORE_HMODEL_HPP
