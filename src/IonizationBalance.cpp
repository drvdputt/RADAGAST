#include "IonizationBalance.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "TemplatedUtils.h"

#include <algorithm>
#include <iostream>

using namespace std;

double Ionization::solveBalance(double nH, double T, const Array& frequencyv,
                                const Array& specificIntensityv)
{
	double integral = photoRateCoeff(frequencyv, specificIntensityv);
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

	// TODO: find better solution for this
	// For very strong radiation fields, some numerical problems can appear... so we cap to 1.
	return min(ne / nH, 1.);
}

double Ionization::photoRateCoeff(const Array& frequencyv, const Array& specificIntensityv)
{
	size_t nFreq = frequencyv.size();
	size_t iThres = TemplatedUtils::index<double>(THRESHOLD, frequencyv);
	Array integrand(nFreq);
	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] *
		               crossSection(frequencyv[n]);
	double integral = Constant::FPI / Constant::PLANCK *
	                  TemplatedUtils::integrate<double>(frequencyv, integrand);
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

		double Fy = (x - 1) * (x - 1) * pow(y, 0.5 * P - 5.5) *
		            pow(1 + sqrt(y / ya), -P);

		return sigma0 * Fy * 1e-18;
	}
}

double Ionization::collisionalRateCoeff(double T)
{
	// 1991-Scholz equations 9 and 10 and table 2
	const vector<double> av{-9.61443e1,  3.79523e1,  -7.96885,   8.83922e-1,
	                        -5.34513e-2, 1.66344e-3, -2.08888e-5};
	double bigGamma = exp(TemplatedUtils::evaluatePolynomial(log(T), av));
	return bigGamma * exp(-Constant::PLANCK * THRESHOLD / Constant::BOLTZMAN / T);
}

double Ionization::recombinationRateCoeff(double T)
{
	double a = 7.982e-11;
	double b = 0.7480;
	double T0 = 3.148;
	double T1 = 7.036e5;

	return a / (sqrt(T / T0) * pow(1 + sqrt(T / T0), 1 - b) * pow(1 + sqrt(T / T1), 1 + b));
}

double Ionization::heating(double np, double ne, double T, const Array& frequencyv,
                           const Array& specificIntensityv)
{
	// Use formula 3.2 from Osterbrock
	double numberOfIonizations = np * ne * recombinationRateCoeff(T);

	size_t nFreq = frequencyv.size();
	Array integrand(nFreq);
	size_t iThres = TemplatedUtils::index<double>(THRESHOLD, frequencyv);

	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] *
		               (frequencyv[n] - THRESHOLD) * crossSection(frequencyv[n]);
	double topIntegral = Constant::PLANCK *
	                     TemplatedUtils::integrate<double>(frequencyv, integrand);

	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] *
		               crossSection(frequencyv[n]);

	// The denominator comes from isolating n_0 from the balance equation, and now also includes
	// the collisional term (top and bottom have been multiplied with h / 4pi, see 3.1, hence
	// the extra
	// factors)
	double bottom = TemplatedUtils::integrate<double>(frequencyv, integrand) +
	                Constant::PLANCK / Constant::FPI * ne * collisionalRateCoeff(T);

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
