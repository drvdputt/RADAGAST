#ifndef _PHOTOELECTRICHEATING_H_
#define _PHOTOELECTRICHEATING_H_

#include "Constants.h"

#include <vector>
#include <cstring>

// TODO: convert everything to frequency units, or write interface

class PhotoelectricHeatingRecipe
{
public:
	void setG0(double G0){_G0 = G0;}
	void setGasTemperature(double T){_gasTemperature = T;}

	double yieldFunctionTest() const;

	// Makes a plot of the heating efficiency in function of the grain size. Saved as a two-column file using the given filename
	double heatingRateTest(std::string filename) const;

	double chargeBalanceTest() const;

private:
	// Functions to calculate the heating rate according to the recipe by Weingartner and Draine (2001)

	double ionizationPotential(double a, int Z) const;

	void chargeBalance(double a, const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
			const std::vector<double>& isrf, int& resultZmax, int& resultZmin,
			std::vector<double>& resultfZ) const;

	double heatingRateAZ(double a, int Z, const std::vector<double>& energyDensity_lambda,
			const std::vector<double>& Qabs, const std::vector<double>& isrf) const;

	double heatingRateA(double a, const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
			const std::vector<double>& energyDensity_lambda) const;

	double heatingRate(const std::vector<double>& wavelengthv, const std::vector<double>& Qabs,
			const std::vector<double>& isrf) const;

	double emissionRate(double a, int Z, const std::vector<double>& wavelengthv,
			const std::vector<double>& Qabs, const std::vector<double>& isrf) const;

	double energyIntegral(double Elow, double Ehigh, double Emin, double Emax) const;

	double yield(double a, int Z, double hnuDiff, double Elow, double Ehigh) const;

	int minimumCharge(double a) const;

	double stickingCoefficient(double a, int Z, int z_i) const;

	double collisionalChargingRate(double a, int Z, int z_i, double m_i) const;

	double lambdaTilde(double tau, double ksi) const;

	double recombinationCoolingRate(double a, const std::vector<double>& fZ, int Zmin) const;

	// Gas properties
	const double _hydrogenDensity{2.5e1};
	const double _ionizationFraction{3.e-4};
	const double _electronDensity{_hydrogenDensity * _ionizationFraction};
	double _gasTemperature{1.e3};

	// Grain properties to use for tests (more detailes in generateQabs)
	const bool _carbonaceous{true};
	const double _workFunction{_carbonaceous ? 4.4 / Constant::ERG_EV : 8 / Constant::ERG_EV};

	// Radiation field to use for test (more details in generateISRF)
	const double _Tc{3.e4};
	double _G0{2.45e0};

	// Wavelength grid to use for tests
	const size_t _nWav{200};
	const double _minWav{0.0912 * Constant::UM_CM}; // cutoff at 13.6 eV
	const double _maxWav{1000 * Constant::UM_CM};

	// Data to use for tests
	void readQabs() const;
	std::vector<double> generateQabs(double a, const std::vector<double>& wavelengthv) const;
};

#endif
