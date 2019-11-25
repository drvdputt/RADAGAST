#include "SimpleH2.hpp"
#include "GasStruct.hpp"
#include "RadiationFieldTools.hpp"

void SimpleH2::solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
                     double h2form)
{
	_nH2 = n;

	// double Rpump = 3.4e-10 * beta(tau) * G_0 * exp(-2.5 * Av); I assume that G_0 *
	// exp(-2.5 * Av) is just the attenuated radiation field. We should calculate G here
	// directly from the given specific intensity.
	_g = RadiationFieldTools::gHabing(specificIntensity);

	// We might needs some value for tau / beta to describe self-shielding. For now, use 1.
	double beta = 1;

	// equation A8 and a factor 1e-1 (from the text: 10% pumps end up in dissociation)
	double Rpump = 3.4e-10 * beta * _g;
	_solomonDissoc = .1 * Rpump;
	// equation A12
	_directDissoc = 1e-11 * _g;

	constexpr double Rdecay = 2e-7; // s-1

	// Equation A13 and A14, for downward collisions --> deexcitation heating.
	double T = gas._t;
	double sqrtT = std::sqrt(T);
	double colH = 1.e-12 * sqrtT * std::exp(-1000. / T) // cm3 s-1
	              * gas._sv.nH();
	double colH2 = 1.4e-12 * sqrtT * std::exp(-18100 / (T + 1200)) // cm3 s-1
	               * gas._sv.nH2();

	// 90% of FUV pumps end up in vib-rot excited states. The latter can decay radiatively
	// or collisionally (eq. A14 + text immediatly below)
	double ratio_ns_ng = .9 * Rpump / (Rdecay + colH + colH2);

	// ns = ratio * ng = ratio * (n - ns) ==> ns = ratio * n / (1 + ratio)
	_nH2s = ratio_ns_ng * _nH2 / (1 + ratio_ns_ng);
	_nH2g = _nH2 - _nH2s;

	// Equation A9, for dissociation heating
	_gamma3 = 1.36e-23 * _nH2 * beta * _g;

	// Energy of the psuedo level. See text below eq. A11 in TH85
	constexpr double Es = 2.6 / Constant::ERG_EV;
	_gamma4 = (colH + colH2) * _nH2s * Es;
}

double SimpleH2::dissociationRate(const Spectrum&) const
{
	// In TH85, only H2* can be dissociated directly
	return _solomonDissoc + _nH2s / _nH2 * _directDissoc;
}

double SimpleH2::dissociationHeating(const Spectrum&) const { return _gamma3; }

double SimpleH2::netHeating() const
{
	// TODO: H2 line cooling per H2 molecule (mostly J-changing collisions). Suggestion:
	// Glover & Abel 2008, MNRAS, 388, 1627; Sections 2.3.1-2.3.6; this is used in cloudy.
	return _gamma4;
}

double SimpleH2::orthoPara() const { return .75; }

Array SimpleH2::emissivityv(const Array& eFrequencyv) const
{
	return Array(eFrequencyv.size());
}

Array SimpleH2::opacityv(const Array& oFrequencyv) const { return Array(oFrequencyv.size()); }
