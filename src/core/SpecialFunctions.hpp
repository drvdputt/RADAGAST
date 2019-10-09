#ifndef _SPECIALFUNCTIONS_
#define _SPECIALFUNCTIONS_

#include "Table.hpp"

#include <functional>

/** Defines a number of mathematical functions, some 'special', others not so special. */
namespace SpecialFunctions
{
/** This function returns the Voigt function \f[ H(a,x) = \frac{a}{\pi} \int_{-\infty}^\infty
    \frac{ e^{-u^2}\,{\rm{d}}u }{ (x-u)^2+a^2 } \f] We use an implementation based on
    the ROOT library. This implementation follows an algorithm developed by Humlicek (1982,
    JQSRT, 21, 437) and Wells (1999, JQSRT, 62, 29), and should be accurate to at least 5
    significant digits. */
double voigt(double a, double u);

/** This function returns the Maxwell-Boltzman velocity distribution, evaluated at a velocity v,
    for a gas of temperature T and particles of mass m. Since this is a probability density per
    velocity unit, the unit of the result is s / cm. */
double maxwellBoltzman(double v, double T, double m);

/** This function evaluates Planck's law at a frequency nu, for a photon gas of temperature T.
    The result is the specific intensity in cgs: erg / s / cm^2 / sr / Hz. */
double planck(double nu, double T);

/** The standard form of a Lorentzian profile, a.k.a. Cauchy distribution. x is the excursion,
    and gamma is the width parameter. For the specific case of a line profile, x should be the
    frequency difference (nu - nu0), and gamma should be the half width the line (total decay
    rate / 4 pi). */
double lorentz(double x, double gamma);

/** For a certain value of the Lorentz profile, return the corresponding (positive) x */
double inverse_lorentz(double l, double gamma);

/** The inverse cumulative distribution a.k.a. percentile function of the lorentz profile */
double lorentz_percentile(double p, double gamma);

double gauss(double x, double sigma);
double inverse_gauss(double g, double sigma);

class LookupTable2D
{
public:
	LookupTable2D(std::function<double(double, double)> ff);
	void setup(const Array& xv, const Array& yv);
	double operator()(double x, double y) const;

private:
	std::function<double(double,double)> _ff;
	Array _xv, _yv;
	double _xMin{0}, _xMax{0}, _yMin{0}, _yMax{0};
	Table<2> _fvv;
};
}

#endif /* _SPECIALFUNCTIONS_ */
