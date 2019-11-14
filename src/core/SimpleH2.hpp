#ifndef CORE_SIMPLEH2_HPP
#define CORE_SIMPLEH2_HPP

#include "H2Model.hpp"

class SimpleH2 : H2Model
{
	void solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
	                   double h2form = 0) override;
	double dissociationRate(const Spectrum& specificIntensity) const override;
	double dissociationHeating(const Spectrum& specificIntensity) const override;
	double netHeating() const override;
	double orthoPara() const override;
	Array emissivityv(const Array& eFrequencyv) const override;
	Array opacityv(const Array& oFrequencyv) const override;
};

#endif // CORE_SIMPLEH2_HPP
