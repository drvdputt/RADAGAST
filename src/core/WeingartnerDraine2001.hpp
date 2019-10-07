#ifndef CORE_WEINGARTNERDRAINE2001_H_
#define CORE_WEINGARTNERDRAINE2001_H_
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
double energyIntegral(double Elow, double Ehigh, double Emin, double Emax);

/** Calculates the photoelectric yield according to WD01 equation 12. */
double yield(double a, int Z, double hnu, bool carbonaceous);

/** Implements WD01 equation 23. */
double autoIonizationThreshold(double a, bool carbonaceous);

/** Implements WD01 equation 24. */
int minimumCharge(double a, double Uait);

/** Combines eq 23 and 24 */
int minimumCharge(double a, bool carbonaceous);

/** The sticking coefficient based on equations 1 and 27-30 */
double stickingCoefficient(double a, int Z, int z_i, bool carbonaceous);

/** Draine & Sutin (1987) equations 3.6-3.10. */
double lambdaTilde(double tau, double ksi);

/** Draine & Sutin (1987) equation 2.4a (with nu replaced by ksi in notation). */
double thetaKsi(double ksi);

} /* namespace WD01 */
#endif /* CORE_WEINGARTNERDRAINE2001_H_ */
