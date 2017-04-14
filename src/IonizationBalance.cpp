#include "IonizationBalance.h"
#include "TemplatedUtils.h"

#include <algorithm>
#include <iostream>

using namespace std;

double Ionization::ionizedFraction(double nH, double T, const Array& frequencyv,
                                   const Array& specificIntensityv)
{
	size_t nFreq = frequencyv.size();
	Array integrand(nFreq);
	size_t iThres = TemplatedUtils::index<double>(THRESHOLD, frequencyv);
	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] * crossSection(frequencyv[n]);

	double C = Constant::FPI / Constant::PLANCK *
	           TemplatedUtils::integrate<double>(frequencyv, integrand) / recombinationRate(T);

	// np^2 = (nH - np) * integral / alpha = (nH - np) * C
	// np^2 + C*np - C*nH = 0
	// np = (-C + sqrt(C^2 + 4 * C*nH)) / 2
	double np = (-C + sqrt(C * C + 4 * C * nH)) / 2.;

	return min(np / nH, 1.);
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

double Ionization::recombinationRate(double T)
{
	double a = 7.982e-11;
	double b = 0.7480;
	double T0 = 3.148;
	double T1 = 7.036e5;

	return a / (sqrt(T / T0) * pow(1 + sqrt(T / T0), 1 - b) * pow(1 + sqrt(T / T1), 1 + b));
}

double Ionization::heating(double nH, double f, double T, const Array& frequencyv,
                           const Array& specificIntensityv)
{
	// Use formula 3.2 from Osterbrock
	double numberOfIonizations = f * f * nH * nH * recombinationRate(T);

	size_t nFreq = frequencyv.size();
	Array integrand(nFreq);
	size_t iThres = TemplatedUtils::index<double>(THRESHOLD, frequencyv);

	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] * (frequencyv[n] - THRESHOLD) *
		               crossSection(frequencyv[n]);
	double topIntegral =
	                Constant::PLANCK * TemplatedUtils::integrate<double>(frequencyv, integrand);

	for (size_t n = iThres; n < nFreq; n++)
		integrand[n] = specificIntensityv[n] / frequencyv[n] * crossSection(frequencyv[n]);
	double bottomIntegral = TemplatedUtils::integrate<double>(frequencyv, integrand);

	double result = numberOfIonizations * topIntegral / bottomIntegral;

	cout << "Ionization heating " << result << endl;

	return result;
}

double Ionization::cooling(double nH, double f, double T)
{
	double kT = Constant::BOLTZMAN * T;
	double T_eV = kT * Constant::ERG_EV;
	double ne = nH * f;

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

	double result = ne * ne * kT * recombinationRate(T) * ft;

	cout << "Ionization cooling " << result << endl;

	return result;
}
