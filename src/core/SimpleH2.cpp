#include "SimpleH2.hpp"
#include "GasStruct.hpp"
#include "RadiationFieldTools.hpp"

void SimpleH2::solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
                     double h2form)
{
	// Energy of the psuedo level. See text below eq. A11 in TH85
	constexpr double Es = 2.6 / Constant::ERG_EV;

	// TODO: properly calculate the correct ratio somehow. Question: are we going to have H2
	// and H2* separately in the chemical network? In that case, the big H2 model also needs
	// a way to convert the level populations in to H2 and H2* densities (shouldn't be hard,
	// just use an energy threshold somewhere).
	_nH2 = n;
	_nH2g = _nH2 / 2;
	_nH2s = _nH2 - _nH2g;
	_g = RadiationFieldTools::gHabing(specificIntensity);

	// escape factor (we might need it to describe self shielding (exintinction of LW line
	// photons which can not be resolved in wavelength space)
	double beta = 1;

	// We might needs some value for tau to describe self-shielding. for now, use 1. Tau
	// depends on the H2 column density and the turbulent Doppler line width (in km s-1).
	// s-1 (equation A8 in TH85), with a factor 0.1 (see cloudy source code

	// double Rpump = 3.4e-10 * beta(tau) * G_0 * exp(-2.5 * Av); I assume that G_0 *
	// exp(-2.5 * Av) is just the attenuated radiation field. We should calculate G here
	// directly from the given specific intensity.

	// equation A8 with a factor 1e-1
	_solomonDissoc = 3.4e-11 * beta * _g;
	// equation A12
	_directDissoc = 1e-11 * _g;

	// Equation A9, for dissociation heating
	_gamma3 = 1.36e-23 * _nH2 * beta * _g;

	// Equation A13 and A14, for downward collisions --> deexcitation heating.
	double T = gas._t;
	double sqrtT = std::sqrt(T);
	double colH = 1.e-12 * sqrtT * std::exp(-1000. / T) // cm3 s-1
	              * gas._sv.nH();
	double colH2 = 1.4e-12 * sqrtT * std::exp(-18100 / (T + 1200)) // cm3 s-1
	               * gas._sv.nH2();
	_gamma4 = (colH + colH2) * _nH2s * Es;
}

double SimpleH2::dissociationRate(const Spectrum&) const
{
	// In TH85, only H2* can be dissociated directly
	return _solomonDissoc + _nH2s / _nH2 * _directDissoc;
}

double SimpleH2::dissociationHeating(const Spectrum&) const { return _gamma3; }

double SimpleH2::netHeating() const { return _gamma4; }

double SimpleH2::orthoPara() const { return .75; }

Array SimpleH2::emissivityv(const Array& eFrequencyv) const
{
	return Array(eFrequencyv.size());
}

Array SimpleH2::opacityv(const Array& oFrequencyv) const { return Array(oFrequencyv.size()); }
