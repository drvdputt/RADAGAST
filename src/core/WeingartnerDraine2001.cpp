#include "WeingartnerDraine2001.hpp"
#include "Constants.hpp"
#include "Error.hpp"
#include "Options.hpp"
#include <gsl/gsl_math.h>
#include <cmath>

using namespace std;

double WD01::atomMass(bool carbonaceous)
{
    // Treat all carbonaceous grains as pure carbon -> use atomic mass. I got the number for
    // silicates from silicate_0m010.opc in cloudy (atomic mass of Si is 28, but theres also some
    // O, Mg and Fe mixed in).
    return carbonaceous ? 12.011 * Constant::AMU_CGS : 24.6 * Constant::AMU_CGS;
}

double WD01::bulkDensity(bool carbonaceous)
{
    // Values taken from SKIRT source code for Draine Graphite and Draine Silicate
    return carbonaceous ? 2.24 : 3.0;
}

double WD01::workFunction(bool carbonaceous)
{
    return carbonaceous ? 4.4 / Constant::ERG_EV : 8 / Constant::ERG_EV;
}

double WD01::eMin(double a, int Z)
{
    if (Options::weingartnerdraine2001_vanHoofEmin)
    {
        // replace by van Hoof (2004) eq 1 (or WD06 eq 3)
        if (Z >= -1)
        {
            return 0;
        }
        else
        {
            double ksi = abs(Z + 1);
            return thetaKsi(ksi) * (1 - 0.3 * pow(a / 10 / Constant::ANG_CM, -0.45) * pow(ksi, -0.26));
        }
    }
    else
    {
        // WD01 eq 7
        double e2_a{Constant::ESQUARE / a};
        double Emin = Z >= 0 ? 0 : -(Z + 1) * e2_a / (1 + pow(27. * Constant::ANG_CM / a, 0.75));
        return Emin;
    }
}

double WD01::ionizationPotential(double a, int Z, bool carbonaceous)
{
    double e2_a = Constant::ESQUARE / a;
    double ip_v = (Z + .5) * e2_a;
    if (Z >= 0)
    {
        // use the same expression for carbonaceous and silicate (WD01 eq 2)
        ip_v += workFunction(carbonaceous) + (Z + 2) * e2_a * 0.3 * Constant::ANG_CM / a;
    }
    // For negatively charged grains, different expressions are used for car and sil
    else if (carbonaceous)
    {
        // WD01 eq 4 (using IP(Z < 0) = EA(Z + 1))
        ip_v += workFunction(carbonaceous) - e2_a * 4.e-8 / (a + 7 * Constant::ANG_CM);
    }
    else  // if silicate
    {
        // WD01 eq 5 (using IP(Z < 0) = EA(Z + 1))
        ip_v += 3 / Constant::ERG_EV;
    }
    return ip_v;
}

double WD01::energyIntegral(double Elow, double Ehigh, double Emin)
{
    // This integrates over a parabolic distribution, with zeros at Elow and Ehigh. The
    // integration range is [Emin, Emax = Ehigh]. Emin has to be given, because we do not
    // want to integrate over electrons that do not escape.

    // I do this to make the code more readable. It looks more like my notes this way.
    double Emax = Ehigh;
    double Ediff = Ehigh - Elow;
    double Ediff3 = Ediff * Ediff * Ediff;

    // Compute integral f(E)E dE analytically, with f(E) defined by WD01 eq 10 f(E) is a parabola,
    // and therefore f(E)E is a third order polynomial. Thus the integral of f(E)E dE is a fourth
    // order polynomial: a/4 (max4 - min4) + b/3 (max3 - min3) + c/2 (max2 - min2).
    double Emax2 = Emax * Emax;
    double Emax3 = Emax2 * Emax;
    double Emax4 = Emax3 * Emax;
    double Emin2 = Emin * Emin;
    double Emin3 = Emin2 * Emin;
    double Emin4 = Emin3 * Emin;
    return 6 / Ediff3
           * (-(Emax4 - Emin4) / 4. + (Ehigh + Elow) * (Emax3 - Emin3) / 3. - Elow * Ehigh * (Emax2 - Emin2) / 2.);
}

double WD01::escapingFraction(int Z, double Elow, double Ehigh)
{
    // Calculate y2 from eq 11
    double Ediff = Ehigh - Elow;
    double Ediff2 = Ediff * Ediff;
    double Ediff3 = Ediff2 * Ediff;
    double y2 = Z >= 0 ? Ediff2 * (Ehigh - 3. * Elow) / Ediff3 : 1;
    return y2;
}

double WD01::yield(double a, int Z, double hnuDiff, double Emin, bool carbonaceous)
{
    double the_y1 = y1(a);
    return yield_cached(a, Z, hnuDiff, Emin, carbonaceous, the_y1);
}

double WD01::yield_cached(double a, int Z, double hnuDiff, double Emin, bool carbonaceous, double y1_cached)
{
    // Below photoelectric threshold
    if (hnuDiff < 0) return 0;

    // Energy distribution parameters (parabola that becomes zero at Elow and Ehigh. See WD01 text
    // between eq 10 and 11).
    double Elow, Ehigh;
    if (Z < 0)
    {
        Elow = Emin;
        Ehigh = hnuDiff + Emin;
    }
    else
    {
        Elow = -(Z + 1) * Constant::ESQUARE / a;
        Ehigh = hnuDiff;
    }

    // Calculate y0 from eq 9, 16, 17
    double thetaOverW = hnuDiff;
    if (Z >= 0) thetaOverW += (Z + 1) * Constant::ESQUARE / a;
    thetaOverW /= workFunction(carbonaceous);

    // Different formulae for carbonaceous and silicate
    double y0;
    if (carbonaceous)
    {
        double thetaOverWtothe5th = gsl_pow_5(thetaOverW);
        y0 = 9e-3 * thetaOverWtothe5th / (1 + 3.7e-2 * thetaOverWtothe5th);
    }
    else
    {
        y0 = 0.5 * thetaOverW / (1. + 5. * thetaOverW);
    }
    double y2 = escapingFraction(Z, Elow, Ehigh);
    return y2 * min(y0 * y1_cached, 1.);
}

double WD01::y1(double a)
{
    // Calculate y1 from grain properties and eq 13, 14
    // double imaginaryRefIndex = 1; // should be wavelength-dependent
    // double la = Constant::LIGHT / nu / 4 / Constant::PI / imaginaryRefIndex;

    // value from 1994-Bakes. WD01 uses the one in the comment above, which is more annoying
    // to calculate.
    constexpr double la = 100 * Constant::ANG_CM;
    constexpr double le = 10 * Constant::ANG_CM;
    double beta = a / la;
    double alpha = beta + a / le;
    double beta2 = beta * beta;
    double alpha2 = alpha * alpha;
    double result =
        beta2 / alpha2 * (alpha2 - 2. * alpha - 2. * expm1(-alpha)) / (beta2 - 2. * beta - 2. * expm1(-beta));
    return result;
}

double WD01::autoIonizationThreshold(double a, bool carbonaceous)
{
    // WD01 eq 23
    double aA = a / Constant::ANG_CM;
    return carbonaceous ? 3.9 + 0.12 * aA + 2. / aA : 2.5 + 0.07 * aA + 8. / aA;
}

int WD01::minimumCharge(double a, double Uait)
{
    // WD01 eq 24
    return floor(-Uait / 14.4 * a / Constant::ANG_CM) + 1;
}

int WD01::minimumCharge(double a, bool carbonaceous)
{
    return minimumCharge(a, autoIonizationThreshold(a, carbonaceous));
}

double WD01::estick_positive(double a)
{
    double le = 10. * Constant::ANG_CM;
    double pElasticScatter = .5;
    return (1 - pElasticScatter) * (-expm1(-a / le));
}

double WD01::estick_negative(double a)
{
    // electron mean free path length in grain
    double le = 10. * Constant::ANG_CM;
    // number of carbon atoms
    double NC = 468 * a * a * a / 1.e-21;
    return 0.5 * (-expm1(-a / le)) / (1 + exp(20 - NC));
}

double WD01::stickingCoefficient_cached(double a, int z, int z_i, bool carbonaceous, double estick_cached_positive,
                                        double estick_cached_negative)
{
    // ions
    if (z_i >= 0) return 1;
    // electrons
    else if (z_i == -1)
    {
        // negative and neutral grains
        if (z <= 0)
        {
            if (z > minimumCharge(a, carbonaceous))
                return estick_cached_negative;
            else
                return 0;
        }
        // positive grains
        else
            return estick_cached_positive;
    }
    // more negative
    else
        return 0;
}

double WD01::stickingCoefficient(double a, int z, int z_i, bool carbonaceous)
{
    double the_estick_positive = estick_positive(a);
    double the_estick_negative = estick_negative(a);
    return stickingCoefficient_cached(a, z, z_i, carbonaceous, the_estick_positive, the_estick_negative);
}

double WD01::lambdaTilde(double tau, double ksi)
{
    /* Found in 1987-Draine-Sutin (writing ksi instead of nu, to avoid confusion with
	   frequency). */
    if (ksi < 0)
    {
        return (2. - ksi / tau) * (1. + 1. / sqrt(tau - ksi));
    }
    else if (ksi > 0)
    {
        return (2. + ksi / tau) * (1. + 1. / sqrt(3. / 2. / tau + 3. * ksi)) * exp(-thetaKsi(ksi) / tau);
    }
    else
    {
        return 2. + 3. / 2. * sqrt(Constant::PI / 2. / tau);
    }
}

double WD01::thetaKsi(double ksi)
{
    // Note that this is an approximation. The exact solution is actually a root of an equation
    if (ksi > 0.)
        return ksi / (1. + 1. / sqrt(ksi));
    else
        return 0;
}

double WD01::sigmaPDT(int Z, double hnuDiff)
{
    constexpr double DeltaE = 3 / Constant::ERG_EV;
    double x = hnuDiff / DeltaE;
    double denom = 1 + x * x / 3;
    // WD01 eq 20, with constant factor moved in front of integral
    double sigma_pdt = 1.2e-17 * abs(Z) * x / denom / denom;
    return sigma_pdt;
}
