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
	double _nH2{0};
	// H2 in pseudo ground state
	double _nH2g{0};
	// H2 in vib-rot excited pseudo level
	double _nH2s{0};
	// Radiation field in habing units
	double _g{0};

	// eq. A8
	double _solomonDissoc{0};
	// eq. A12
	double _directDissoc{0};
	// Heating due to solomon dissociation, eq. A9
	double _gamma3;
	// Collisional de-excitation heating, eq. A13
	double _gamma4;
};

#endif // CORE_SIMPLEH2_HPP
