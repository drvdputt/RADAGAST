#include "Ionization.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "TemplatedUtils.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

namespace RADAGAST
{
    double Ionization::solveBalance(double nH, double T, const Spectrum& meanIntensity)
    {
        double integral = photoRateCoeff(meanIntensity);
        double gamma = collisionalRateCoeff(T);
        double alpha = recombinationRateCoeff(T);

        // Solve quadratic equation
        double c2 = 1;
        double c1 = (-nH * gamma + integral) / (alpha + gamma);
        double c0 = -nH * integral / (alpha + gamma);
        double ne = (-c1 + sqrt(c1 * c1 - 4 * c2 * c0)) / 2.;
        // old (radiative only):
        //		double C = Constant::FPI / Constant::PLANCK *
        //                     TemplatedUtils::integrate<double>(frequencyv, integrand) /
        //                     recombinationRateCoeff(T);
        //
        //		// np^2 = (nH - np) * integral / alpha = (nH - np) * C
        //		// np^2 + C*np - C*nH = 0
        //		// np = (-C + sqrt(C^2 + 4 * C*nH)) / 2
        //		double np = (-C + sqrt(C * C + 4 * C * nH)) / 2.;
        //		ne = np;

        //
        DEBUG("Ionization rate is " << ne * recombinationRateCoeff(T) << " s-1" << endl);

        // For very strong radiation fields, some numerical problems can appear... so we cap to 1.
        return min(ne / nH, 1.);
    }

    double Ionization::photoRateCoeff(const Spectrum& meanIntensity)
    {
        const Array& Iv = meanIntensity.valuev();
        const Array& nuv = meanIntensity.frequencyv();

        // start at largest frequency, and stop when we go below threshold
        int i = nuv.size() - 1;

        // when ionizing radiation is not included in wavelength grid, return early
        if (nuv[i] < THRESHOLD) return 0.;

        Array integrandv(nuv.size());
        while (i >= 0 && nuv[i] > THRESHOLD)
        {
            // Start from I_nu -> erg s-1 cm-2 hz-1 sr-1. Then, 4pi I_nu / hnu = s-1 cm-2 hz-1.
            // Then multiply with crossSection -> s-1 hz-1.
            integrandv[i] = Iv[i] / nuv[i] * crossSection(nuv[i]);
            i--;
        }
        // i is now 1 below integration bound

        // integrate over frequency -> s-1
        double rate = Constant::FPI / Constant::PLANCK
                      * TemplatedUtils::integrate<double>(nuv, integrandv, i + 1, nuv.size() - 1);
        return rate;
    }

    double Ionization::crossSection(double frequency)
    {
        if (frequency < THRESHOLD)
        {
            return 0;
        }
        else
        {
            double E0 = 4.298e-1;
            double sigma0 = 5.475e4;
            double ya = 3.288e1;
            double P = 2.963;

            double x = Constant::PLANCK * frequency / (E0 * Constant::EV);
            double y = x;

            double Fy = (x - 1) * (x - 1) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);

            return sigma0 * Fy * 1e-18;
        }
    }

    double Ionization::collisionalRateCoeff(double T)
    {
        // 1991-Scholz equations 9 and 10 and table 2
        const vector<double> av{-9.61443e1, 3.79523e1, -7.96885, 8.83922e-1, -5.34513e-2, 1.66344e-3, -2.08888e-5};
        double bigGamma = exp(TemplatedUtils::evaluatePolynomial(log(T), av));
        return bigGamma * exp(-Constant::PLANCK * THRESHOLD / Constant::BOLTZMAN / T);
    }

    double Ionization::recombinationRateCoeff(double T)
    {
        // These coefficients are from Verner and Ferland (1996). Those from Badnell (2006) are
        // slightly different.
        double a = 7.982e-11;
        double b = 0.7480;
        double T0 = 3.148;
        double T1 = 7.036e5;

        return a / (sqrt(T / T0) * pow(1 + sqrt(T / T0), 1 - b) * pow(1 + sqrt(T / T1), 1 + b));
    }

    double Ionization::heatingPerH(const Spectrum& meanIntensity)
    {
        const Array& Iv = meanIntensity.valuev();
        const Array& nuv = meanIntensity.frequencyv();

        // stop when we go below threshold
        int i = nuv.size() - 1;
        // when ionizing radiation is not included in wavelength grid, return early
        if (nuv[i] < THRESHOLD) return 0.;

        Array integrandv(nuv.size());
        while (i >= 0 && nuv[i] > THRESHOLD)
        {
            // Start from I_nu -> erg s-1 cm-2 hz-1 sr-1. Then, h(nu - threshold) * 4pi I_nu /
            // hnu = (nu - threshold) * 4pi I_nu / nu -> erg cm-2 hz-1. Then multiply with
            // crossSection -> erg s-1 hz-1.
            integrandv[i] = (nuv[i] - THRESHOLD) * Constant::FPI * Iv[i] / nuv[i] * crossSection(nuv[i]);
            i--;
        }
        // i is now 1 below integration bound

        // integrate over frequency -> erg s-1
        double heatPerH = TemplatedUtils::integrate<double>(nuv, integrandv, i + 1, nuv.size() - 1);
        return heatPerH;
    }

    double Ionization::cooling(double nH, double np, double ne, double T)
    {
        // Kinetic energy lost through radiative recombination
        double kT = Constant::BOLTZMAN * T;
        double T_eV = kT / Constant::EV;

        // use fit from 2017-Mao ft = beta / alpha cooling is rate coefficient is kT * beta = kT *
        // alpha * ft

        double a0 = 8.655e0;
        double b0 = 5.432e-1;
        //	double c0 = 0;
        double a1 = 1.018e1;
        double b1 = 5.342e-1;
        //	double a2 = 0;
        //	double b2 = 0;
        //	double ft = a0 * pow(T_eV, -b0 - c0 * log10(T_eV)) * (1 + a2 * pow(T_eV, -b2)) /
        //              (1 + a1 * pow(T_eV, -b1));
        double ft = a0 * pow(T_eV, -b0) / (1 + a1 * pow(T_eV, -b1));

        double result = np * ne * kT * recombinationRateCoeff(T) * ft;

        // Kinetic energy lost through collisional ionization (rate * ionization potential)
        result += Constant::PLANCK * Ionization::THRESHOLD * ne * nH * collisionalRateCoeff(T);

        return result;
    }
}
