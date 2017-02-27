#ifndef _PHOTOELECTRICHEATING_H_
#define _PHOTOELECTRICHEATING_H_

#include <vector>
#include <cstring>

// TODO: convert everything to frequency units, or write interface

namespace Photoelectric
{
// Functions to calculate the heating rate according to the recipe by Weingartner and Draine (2001)

double ionizationPotential(double a, int Z);

double heatingRateAZ(double a, int Z, const std::vector<double>& energyDensity_lambda,
		const std::vector<double>& Qabs, const std::vector<double>& isrf);

double heatingRateA(double a, const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
		const std::vector<double>& energyDensity_lambda);

double heatingRate(const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
		const std::vector<double>& isrf);

double emissionRate(double a, int Z, const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
		const std::vector<double>& isrf);

double energyIntegral(double Elow, double Ehigh, double Emin, double Emax);

double yield(double a, int Z, double hnuDiff, double Elow, double Ehigh);

int minimumCharge(double a);

double stickingCoefficient(double a, int Z, int z_i);

double collisionalChargingRate(double a, int Z, int z_i, double m_i);

double lambdaTilde(double tau, double ksi);

double recombinationCoolingRate(double a, const std::vector<double>& fZ, int Zmin);

void chargeBalance(double a, const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
		const std::vector<double>& isrf, int& resultZmax, int& resultZmin,
		std::vector<double>& resultfZ);

double yieldFunctionTest();

double heatingRateTest(std::string filename);

double chargeBalanceTest();

}

#endif
