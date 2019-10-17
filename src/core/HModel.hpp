#ifndef CORE_HMODEL_HPP
#define CORE_HMODEL_HPP

#include "Array.hpp"
#include "LevelSolution.hpp"

class HData;
class GasStruct;

class HModel
{
public:
	/** Pass a pointer to the HData object at construction, so that the constant H data and
	    functions can be accessed. */
	HModel(const HData* hData) : _hData{hData} {}

	/** Solve the H levels, and store them in this object. */
	void solve(const GasStruct& gas, const Spectrum& specificIntensity);

	/** This function returns the line emission spectrum + the continuum emitted by the
	    2s-1s two-photon process. */
	Array emissivity(const Array eFrequencyv) const;

	/** This function returns the line opacity. TODO: maybe include bound-free cross section
	    that depends on the level populations. */
	Array opacityv(const Array oFrequencyv) const;

	/** From the level populations, calculate the net heating-cooling balance by
	    (de-)excitation */
	double netHeating() const;

private:
	const HData* _hData;
	LevelSolution _levelSolution;
};

#endif // CORE_HMODEL_HPP
