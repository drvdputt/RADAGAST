#include "TemplatedUtils.h"

#include "IonizationBalance.h"

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

	return np / nH;
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
