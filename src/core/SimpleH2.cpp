#include "SimpleH2.hpp"

void SimpleH2::solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
                     double h2form)
{
}

double SimpleH2::dissociationRate(const Spectrum& specificIntensity) const
{
	// double Rpump = 3.4e-10 * beta(tau) * G_0 * exp(-2.5 * Av); I assume that G_0 *
	// exp(-2.5 * Av) is just the attenuated radiation field. We should calculate G here
	// directly from the given specific intensity.

	// We might needs some value for tau to describe self-shielding. for now, use 1. Tau
	// depends on the H2 column density and the turbulent Doppler line width (in km s-1).
	double G = 0;
	double beta = 1;
	double Rpump = 3.4e-10 * beta * G; // s-1 (equation A8 in TH85)
	return 0;
}

double SimpleH2::dissociationHeating(const Spectrum& specificIntensity) const
{
	double G = 0;
	double beta = 1;
	double Gamma3 = 1.36e-23 * nH2 * beta * G;
	return Gamma3;
}

double SimpleH2::netHeating() const { return 0; }

double SimpleH2::orthoPara() const { return 0; }

Array SimpleH2::emissivityv(const Array& eFrequencyv) const
{
	return Array(eFrequencyv.size());
}

Array SimpleH2::opacityv(const Array& oFrequencyv) const { return Array(oFrequencyv.size()); }
