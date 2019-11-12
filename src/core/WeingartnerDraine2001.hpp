#ifndef CORE_WEINGARTNERDRAINE2001_HPP
#define CORE_WEINGARTNERDRAINE2001_HPP

namespace WD01
{
/* Functions to calculate the heating rate according to the recipe by Weingartner and Draine
   (2001) */

/** Return the work function as mentioned in the paper, for carbonaceous/graphitic grains
    (carbonaceous == true) or silicate grains (carbonaceous == false). */
double workFunction(bool carbonaceous);

/** Formula for Emin (WD01 eq 7 replaced by van Hoof (2004) eq 1). */
double eMin(double a, int Z);

/** Implements WD01 equations 2, 4 and 5. */
double ionizationPotential(double a, int Z, bool carbonaceous);

/** Calculates the integral over the energy in WD01 equation 39. */
double energyIntegral(double Elow, double Ehigh, double Emin);

/** Calculate y2 from eq 11, the normalization of the averag energy integral. */
double escapingFraction(int Z, double Elow, double Ehigh);

/** Calculates the photoelectric yield according to WD01 equation 12. y1 is expensive to
    calculate (needs two expm1 evaluations), and is passed as an argument as an optimization. */
double yield_cached(double a, int Z, double hnuDiff, double Emin, bool carbonaceous,
                    double y1_cached);

/** The unoptimized version of the yield function (useful for plotting). */
double yield(double a, int Z, double hnuDiff, double Emin, bool carbonaceous);

/** Implements equations 13 and 14 */
double y1(double a);

/** Implements WD01 equation 23. */
double autoIonizationThreshold(double a, bool carbonaceous);

/** Implements WD01 equation 24. */
int minimumCharge(double a, double Uait);

/** Combines eq 23 and 24 */
int minimumCharge(double a, bool carbonaceous);

/** The sticking coefficient based on equations 1 and 27-30 */
double estick_positive(double a);
double estick_negative(double a);
double stickingCoefficient_cached(double a, int z, int z_i, bool carbonaceous,
                                  double estick_cached_positive, double estick_cached_negative);
double stickingCoefficient(double a, int z, int z_i, bool carbonaceous);

/** Draine & Sutin (1987) equations 3.6-3.10. */
double lambdaTilde(double tau, double ksi);

/** Draine & Sutin (1987) equation 2.4a (with nu replaced by ksi in notation). */
double thetaKsi(double ksi);

/** Photodetachment cross section (equation 20) */
double sigmaPDT(int Z, double hnuDiff);

} /* namespace WD01 */

#endif // CORE_WEINGARTNERDRAINE2001_HPP
