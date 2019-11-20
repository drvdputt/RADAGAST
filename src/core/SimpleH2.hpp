#ifndef CORE_SIMPLEH2_HPP
#define CORE_SIMPLEH2_HPP

#include "H2Model.hpp"

class SimpleH2 : H2Model
{
public:
	void solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
	           double h2form = 0) override;
	double dissociationRate(const Spectrum& specificIntensity) const override;
	double dissociationHeating(const Spectrum& specificIntensity) const override;
	double netHeating() const override;
	double orthoPara() const override;
	Array emissivityv(const Array& eFrequencyv) const override;
	Array opacityv(const Array& oFrequencyv) const override;

private:
	// Try using Burton (1990) appendix A. It has some differences from the Tielens (1985)
	// model, mainly for the H2 line emission and H2 vibrational heating models. I assume
	// the dissociation was still the same as a in Tielens (1985).

	// There are four (v = 0, 1, 2, 3) vibrational levels, and v = 6 is treated as a
	// pseudolevel. The density of the latter is denoted as n^*_{H2}
	double _nH2Star;
	const simpleH2Data* _simpleH2Data;
	LevelSolution _levelSolution;
};

#endif // CORE_SIMPLEH2_HPP
