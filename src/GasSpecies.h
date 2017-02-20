#ifndef _GASSPECIES_H_
#define _GASSPECIES_H_

#include "TwoLevel.h"

#include <vector>

class GasSpecies {
public:
	GasSpecies(const std::vector<double>& wavelengthv);

	void solveBalance(double n, const std::vector<double>& isrf);

	std::vector<double> emissivity() const;
	std::vector<double> opacity() const;

private:
	// Calculates the ionization fraction and level populations for a certain electron temperature and isrf
	void calculateDensities(double T, const std::vector<double>& isrf);

	// Calculate the emissivity, interpolated from data
	std::vector<double> continuumEmissivity(double T);

	// To be set in constructor
	const std::vector<double>& _wavelengthv;

	// Data to be loaded in constructor body
	// First index is for frequency, second for temperature
	std::vector<std::vector<double>> _gammaDaggervv;
	// The accompanying temperature grid
	std::vector<double> _logTemperaturev;

	// To be set on invocation of solveBalance()
	double _n;

	// Results of solveBalance()
	double _ionizedFraction;
	TwoLevel _levels;
};

#endif /* _GASSPECIES_H_ */
