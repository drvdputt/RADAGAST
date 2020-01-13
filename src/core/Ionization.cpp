#include "Ionization.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "TemplatedUtils.hpp"
#include <algorithm>
#include <iostream>

using namespace std;

double Ionization::solveBalance(double nH, double T, const Spectrum& specificIntensity)
{
    double integral = photoRateCoeff(specificIntensity);
    double gamma = collisionalRateCoeff(T);
    double alpha = recombinationRateCoeff(T);

    // Solve quadratic equation
    double c2 = 1;
    double c1 = (-nH * gamma + integral) / (alpha + gamma);
    double c0 = -nH * integral / (alpha + gamma);
    double ne = (-c1 + sqrt(c1 * c1 - 4 * c2 * c0)) / 2.;
    // old (radiative only):
    //		double C = Constant::FPI / Constant::PLANCK *
    //		           TemplatedUtils::integrate<double>(frequencyv, integrand) /
    //		           recombinationRateCoeff(T);
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

double Ionization::photoRateCoeff(const Spectrum& specificIntensity)
{
    // Just integrate over the points of the input spectrum here, since the photoionization
    // cross section is smooth
    const Array& nuv = specificIntensity.frequencyv();
    const Array& vv = specificIntensity.valuev();

    // Integrate over the points [THRESHOLD, nuv[iThres], nuv[iThres + 1], ...]. Including
    // the threshold in the integration is a good correction if the grid is coarse, but
    // ultimately, the user needs to make sure that the grid is fine enough near the peak/edge of
    // the cross section.
    size_t iThres = TemplatedUtils::index<double>(THRESHOLD, nuv);
    size_t numPoints = nuv.size() - iThres + 1;
    Array xv(numPoints);
    Array integrandv(numPoints);

    xv[0] = THRESHOLD;
    copy(begin(nuv) + iThres, end(nuv), begin(xv) + 1);

    // radiation field / nu * cross section

    // interpolate for this point
    integrandv[0] = specificIntensity.evaluate(xv[0]) / xv[0] * crossSection(xv[0]);
    for (size_t i = 1; i < numPoints; i++)
        // use the raw spectrum data for the rest
        integrandv[i] = vv[iThres + i - 1] / xv[i] * crossSection(xv[i]);

    double integral = Constant::FPI / Constant::PLANCK * TemplatedUtils::integrate<double>(xv, integrandv);
    return integral;
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

        double x = Constant::PLANCK * frequency * Constant::ERG_EV / E0;
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

double Ionization::heating(double np, double ne, double T, const Spectrum& specificIntensity)
{
    // Use formula 3.2 from Osterbrock
    double numberOfIonizations = np * ne * recombinationRateCoeff(T);

    auto nuv = specificIntensity.frequencyv();
    auto iv = specificIntensity.valuev();

    size_t nFreq = nuv.size();
    Array integrand(nFreq);
    size_t iThres = TemplatedUtils::index<double>(THRESHOLD, nuv);

    // Integrand I_nu / nu * sigma
    for (size_t i = iThres; i < nFreq; i++) integrand[i] = iv[i] / nuv[i] * crossSection(nuv[i]);

    // The denominator comes from isolating n_0 from the balance equation, and now also
    // includes the collisional term (top and bottom have been multiplied with h / 4pi, see
    // 3.1, hence the extra factors)
    double bottom = TemplatedUtils::integrate<double>(nuv, integrand)
                    + Constant::PLANCK / Constant::FPI * ne * collisionalRateCoeff(T);

    // Now multiply the original integrand with (nu - nu_0). Top / bottom will give the
    // average energy left over after ionization.
    for (size_t i = iThres; i < nFreq; i++) integrand[i] *= (nuv[i] - THRESHOLD);

    double topIntegral = Constant::PLANCK * TemplatedUtils::integrate<double>(nuv, integrand);

    // The heating is hence the total ionization rate (equal to the recombination rate since
    // we assume equilibrium), times the average energy of the release electron.
    double result = numberOfIonizations * topIntegral / bottom;

    DEBUG("Ionization heating " << result << endl);

    return result;
}

double Ionization::cooling(double nH, double np, double ne, double T)
{
    // Kinetic energy lost through radiative recombination
    double kT = Constant::BOLTZMAN * T;
    double T_eV = kT * Constant::ERG_EV;

    // use fit from 2017-Mao
    // ft = beta / alpha
    // cooling is rate coefficient is kT * beta = kT * alpha * ft

    double a0 = 8.655e0;
    double b0 = 5.432e-1;
    //	double c0 = 0;
    double a1 = 1.018e1;
    double b1 = 5.342e-1;
    //	double a2 = 0;
    //	double b2 = 0;
    //	double ft = a0 * pow(T_eV, -b0 - c0 * log10(T_eV)) * (1 + a2 * pow(T_eV, -b2)) /
    //	            (1 + a1 * pow(T_eV, -b1));
    double ft = a0 * pow(T_eV, -b0) / (1 + a1 * pow(T_eV, -b1));

    double result = np * ne * kT * recombinationRateCoeff(T) * ft;

    // Kinetic energy lost through collisional ionization (rate * ionization potential)
    result += Constant::PLANCK * Ionization::THRESHOLD * ne * nH * collisionalRateCoeff(T);

    return result;
}
