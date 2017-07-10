#ifndef _SPECIALFUNCTIONS_
#define _SPECIALFUNCTIONS_

/* Defines a number of mathematical functions, some 'special', others not so special. */
namespace SpecialFunctions
{
/** This function returns the Voigt function \f[ H(a,x) = \frac{a}{\pi} \int_{-\infty}^\infty \frac{
    {\text{e}}^{-u^2}\,{\text{d}}u }{ (x-u)^2+a^2 } \f] We use an implementation based on the ROOT
    library. This implementation follows an algorithm developed by Humlicek (1982, JQSRT, 21, 437)
    and Wells (1999, JQSRT, 62, 29), and should be accurate to at least 5 significant digits. */
double voigt(double a, double u);

/** This function returns the Maxwell-Boltzman velocity distribution, evaluated at a velocity v, for
    a gas of temperature T and particles of mass m. Since this is a probability density per velocity
    unit, the unit of the result is s / cm. */
double maxwellBoltzman(double v, double T, double m);

/** This function evaluates Planck's law at a frequency nu, for a photon gas of temperature T. The
    result is the specific intensity in cgs: erg / s / cm^2 / sr / Hz. */
double planck(double nu, double T);
}

#endif /* _SPECIALFUNCTIONS_ */
