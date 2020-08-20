#ifndef _SPECIALFUNCTIONS_
#define _SPECIALFUNCTIONS_

#include "Table.hpp"
#include <functional>

namespace RADAGAST
{
    /** Defines a number of simple mathematical and physical functions, not bound to a specific
        piece of implementation. */
    namespace Functions
    {
        /** This function returns the Voigt function \f[ H(a,x) = \frac{a}{\pi}
            \int_{-\infty}^\infty \frac{ e^{-u^2}\,{\rm{d}}u }{ (x-u)^2+a^2 } \f] We use an
            implementation based on the ROOT library. This implementation follows an algorithm
            developed by Humlicek (1982, JQSRT, 21, 437) and Wells (1999, JQSRT, 62, 29), and
            should be accurate to at least 5 significant digits. */
        double voigt(double a, double u);

        double pseudoVoigt(double x, double sigma_gauss, double gamma_lorentz);

        /** This function returns the Maxwell-Boltzman velocity distribution, evaluated at a
            velocity v, for a gas of temperature T and particles of mass m. Since this is a
            probability density per velocity unit, the unit of the result is s / cm. */
        double maxwellBoltzman(double v, double T, double m);

        /** Mean particle velocity of Maxwell distribution. Typical use is a v-factor for the
            rate of a collisional process (e.g. grain charging, grain H2 formation). Note that
            this is not the same as the RMS velocity. */
        double meanThermalVelocity(double T, double m);

        /** Width of a Gaussian thermal velocity distribution. Note the difference with the mean
            thermal velocity above: the latter is the mean of the absolute value of the velocity
            (abs(v) is Maxwell-distributed), while this function gives the width of the v
            distribution, which is Gaussian. */
        double thermalVelocityWidth(double T, double m);

        /** This function evaluates Planck's law at a frequency nu, for a photon gas of temperature
            T. The result is the specific intensity in cgs: erg / s / cm^2 / sr / Hz. */
        double planck(double nu, double T);

        /** The standard form of a Lorentzian profile, a.k.a. Cauchy distribution. x is the
            excursion, and gamma is the width parameter. For the specific case of a line profile, x
            should be the frequency difference (nu - nu0), and gamma should be the half width the
            line (total decay rate / 4 pi). */
        double lorentz(double x, double gamma);

        /** For a certain value of the Lorentz profile, return the corresponding (positive) x */
        double inverse_lorentz(double l, double gamma);

        /** The inverse cumulative distribution a.k.a. percentile function of the lorentz
            profile */
        double lorentz_percentile(double p, double gamma);

        double gauss(double x, double sigma);
        double inverse_gauss(double g, double sigma);
    }
}
#endif /* _SPECIALFUNCTIONS_ */
