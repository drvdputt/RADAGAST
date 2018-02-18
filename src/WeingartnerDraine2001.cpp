#include "WeingartnerDraine2001.h"
#include "Constants.h"
#include "Error.h"

#include <cmath>

using namespace std;

double WD01::workFunction(bool carbonaceous)
{
	return carbonaceous ? 4.4 / Constant::ERG_EV : 8 / Constant::ERG_EV;
}

double WD01::eMin(double a, int Z)
{
// TODO: use this as a default instead. Before I can do that, the charge calculation algorithm
// needs to be made more robust.
// #define VANHOOF_EMIN
#ifdef VANHOOF_EMIN
	// replace by van Hoof (2004) eq 1
	if (Z >= -1)
		return 0;
	else
	{
		double ksi{-static_cast<double>(Z) - 1};
		return thetaKsi(ksi) *
		       (1 - 0.3 * pow(a / 10 / Constant::ANG_CM, -0.45) * pow(ksi, -0.26));
	}
#else
	// WD01 eq 7
	double e2_a{Constant::ESQUARE / a};
	double Emin = Z >= 0 ? 0
	                     : -(Z + 1) * e2_a / (1 + pow(27. * Constant::ANG_CM / a, 0.75));
	return Emin;
#endif
}

double WD01::ionizationPotential(double a, int Z, bool carbonaceous)
{
	double e2_a = Constant::ESQUARE / a;
	double ip_v = (Z + .5) * e2_a;
	if (Z >= 0)
	{
		// use the same expression for carbonaceous and silicate
		ip_v += workFunction(carbonaceous) +
		        (Z + 2) * e2_a * 0.3 * Constant::ANG_CM / a; // WD01 eq 2
	}
	// For negatively charged grains, different expressions are used for car and sil
	else if (carbonaceous)
	{
		ip_v += workFunction(carbonaceous) -
		        e2_a * 4.e-8 / (a + 7 * Constant::ANG_CM); // WD01 eq 4
	}
	else // if silicate
	{
		ip_v += 3 / Constant::ERG_EV; // WD01 eq 5
	}
	return ip_v;
}

double WD01::energyIntegral(double Elow, double Ehigh, double Emin, double Emax)
{
	double Ediff = Ehigh - Elow;
	double Ediff3 = Ediff * Ediff * Ediff;

	/* Compute integral f(E)E dE analytically, with f(E) defined by WD01 eq 10 f(E) is a
	   parabola, and therefore f(E)E is a third order polynomial.  Thus the integral of f(E)E dE
	   is a fourth order polynomial: a/4 (max4 - min4) + b/3 (max3 - min3) + c/2 (max2
	   -min2). */
	double Emax2 = Emax * Emax;
	double Emin2 = Emin * Emin;
	return 6 / Ediff3 *
	       (-(Emax2 * Emax2 - Emin2 * Emin2) / 4. +
	        (Ehigh + Elow) * (Emax2 * Emax - Emin2 * Emin) / 3. -
	        Elow * Ehigh * (Emax2 - Emin2) / 2.);
}

double WD01::yield(double a, int Z, double hnu, bool carbonaceous)
{
	double ip_v = ionizationPotential(a, Z, carbonaceous);
	double Emin = eMin(a, Z);
	double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin;
	double hnuDiff = hnu - hnu_pet;
	// Below photoelectric threshold
	if (hnuDiff < 0)
		return 0;

	double Emax = hnuDiff + Emin;

	// Compute yield (y2, y1, y0, and finally Y)

	// WD01 text between eq 10 and 11
	double Elow = Z < 0 ? Emin : -(Z + 1) * Constant::ESQUARE / a;
	double Ehigh = Z < 0 ? Emax : hnuDiff;
	// Calculate y2 from eq 11
	double Ediff = Ehigh - Elow;
	double y2 = Z >= 0 ? Ehigh * Ehigh * (Ehigh - 3. * Elow) / Ediff / Ediff / Ediff : 1;

	// Calculate y1 from grain properties and eq 13, 14
	// double imaginaryRefIndex = 1; // should be wavelength-dependent
	// double la = Constant::LIGHT / nu / 4 / Constant::PI / imaginaryRefIndex;
	const double la = 100 *
	                  Constant::ANG_CM; // value from 1994-Bakes. WD01 uses the above one
	double beta = a / la;
	const double le = 10 * Constant::ANG_CM;
	double alpha = beta + a / le;

	double alpha2 = alpha * alpha;
	double beta2 = beta * beta;
	double y1 = beta2 / alpha2 * (alpha2 - 2. * alpha - 2. * expm1(-alpha)) /
	            (beta2 - 2. * beta - 2. * expm1(-beta));

	// Calculate y0 from eq 9, 16, 17
	double thetaOverW = hnuDiff;
	if (Z >= 0)
		thetaOverW += (Z + 1) * Constant::ESQUARE / a;
	thetaOverW /= workFunction(carbonaceous);
	double thetaOverWtothe5th = thetaOverW * thetaOverW;
	thetaOverWtothe5th *= thetaOverWtothe5th * thetaOverW;

	double y0;
	// Different formulae for carbonaceous and silicate
	if (carbonaceous)
		y0 = 9e-3 * thetaOverWtothe5th / (1 + 3.7e-2 * thetaOverWtothe5th);
	else
		y0 = 0.5 * thetaOverW / (1. + 5. * thetaOverW);

	return y2 * min(y0 * y1, 1.);
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

double WD01::stickingCoefficient(double a, int Z, int z_i, bool carbonaceous)
{
	// ions
	if (z_i >= 0)
		return 1;
	// electrons
	else if (z_i == -1)
	{
		// electron mean free path length in grain
		double le = 10. * Constant::ANG_CM;
		double pElasticScatter = .5;
		// negative and neutral grains
		if (Z <= 0)
		{
			if (Z > minimumCharge(a, carbonaceous))
			{
				// number of carbon atoms
				double NC = 468 * a * a * a / 1.e-21;
				return 0.5 * (-expm1(-a / le)) / (1 + exp(20 - NC));
			}
			else
				return 0;
		}
		// positive grains
		else
			return (1 - pElasticScatter) * (-expm1(-a / le));
	}
	// more negative
	else
		return 0;
}

double WD01::lambdaTilde(double tau, double ksi)
{
	/* Found in 1987-Draine-Sutin (writing ksi instead of nu, to avoid confuction with
	   frequency). */
	if (ksi < 0)
	{
		return (2. - ksi / tau) * (1. + 1. / sqrt(tau - ksi));
	}
	else if (ksi > 0)
	{
		return (2. + ksi / tau) * (1. + 1. / sqrt(3. / 2. / tau + 3. * ksi)) *
		       exp(-thetaKsi(ksi) / tau);
	}
	else
	{
		return 2. + 3. / 2. * sqrt(Constant::PI / 2. / tau);
	}
}

double WD01::thetaKsi(double ksi)
{
	// Note that this is an approximation. The exact soluation is actually a root of an equation
	return ksi > 0 ? ksi / (1. + 1. / sqrt(ksi)) : 0;
}
