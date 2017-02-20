#ifndef _PHOTOELECTRICHEATING_H_
#define _PHOTOELECTRICHEATING_H_

#include <vector>
#include <cstring>

namespace Photoelectric {
// Functions to calculate the heating rate according to the recipe by Weingartner and Draine (2001)

double ionizationPotential(double a, int Z);

double heatingRateAZ(double a, int Z, std::vector<double>& wavelengthv, std::vector<double>& Qabs, std::vector<double>& isrf);

double heatingRateA(double a, std::vector<double>& wavelengthv, std::vector<double>& Qabs, std::vector<double>& isrf);

double heatingRate(std::vector<double>& wavelengthv, std::vector<double>& Qabs, std::vector<double>& isrf);

double emissionRate(double a, int Z, std::vector<double>& wavelengthv, std::vector<double>& Qabs, std::vector<double>& isrf);

double energyIntegral(double Elow, double Ehigh, double Emin, double Emax);

double yield(double a, int Z, double hnuDiff, double Elow, double Ehigh);

int minimumCharge(double a);

double stickingCoefficient(double a, int Z, int z_i);

double collisionalChargingRate(double a, int Z, int z_i, double m_i);

double lambdaTilde(double tau, double ksi);

double recombinationCoolingRate(double a, const std::vector<double>& fZ, int Zmin);

void chargeBalance(double a, std::vector<double>& wavelengthv, std::vector<double>& Qabs, std::vector<double>& isrf,
                   int& resultZmax, int& resultZmin, std::vector<double>& resultfZ);

double yieldFunctionTest();

double heatingRateTest(std::string filename);

double chargeBalanceTest();

}

#endif
