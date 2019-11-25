#ifndef CORE_SIMPLEH2_HPP
#define CORE_SIMPLEH2_HPP

#include "H2Model.hpp"

class LookupTable;

class SimpleH2 : public H2Model
{
public:
	SimpleH2(const LookupTable* lteCool);
	void solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
	           double h2form = 0) override;
	double dissociationRate(const Spectrum& specificIntensity) const override;
	double dissociationHeating(const Spectrum& specificIntensity) const override;
	double netHeating() const override;
	double orthoPara() const override;
	Array emissivityv(const Array& eFrequencyv) const override;
	Array opacityv(const Array& oFrequencyv) const override;

private:
	double gloverAbel08Cooling(const GasStruct& gas) const;

	// This implements the H2 model of Tielens and Hollenback (1985).

	// Question: are we going to have H2 and H2* separately in the chemical network? In that
	// case, the big H2 model also needs a way to convert the level populations in to H2 and
	// H2* densities (shouldn't be hard, just use an energy threshold somewhere).
	double _nH2{0};
	// H2 in pseudo ground state
	double _nH2g{0};
	// H2 in vib-rot excited pseudo level
	double _nH2s{0};
	// ortho fraction
	double _ortho{.75};
	// para fraction
	double _para{1. - _ortho};
	// Radiation field in habing units
	double _g{0};

	// eq. A8
	double _solomonDissoc{0};
	// eq. A12
	double _directDissoc{0};
	// Heating due to solomon dissociation, eq. A9
	double _gamma3{0};
	// Collisional de-excitation heating, eq. A13, due to UV pumping (first to another E
	// level, then back to X, with higher J or v)
	double _gamma4{0};
	// Cooling due to upward excitation by H, H2, e-, p+.
	double _collisionalExcitationCooling{0};
	// Cooling in the high density limit (LTE, precalculated).
	const LookupTable* _lteCool;
};

#endif // CORE_SIMPLEH2_HPP
