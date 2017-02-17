#ifndef _GASSPECIES_H_
#define _GASSPECIES_H_

#include "TwoLevel.h"

#include <vector>

class GasSpecies {
public:
	GasSpecies(double n, const std::vector<double>& isrf, const std::vector<double>& wavelength);

	void solveBalance();

	std::vector<double> emissivity() const;
	std::vector<double> opacity() const;

private:
	// Calculates the ionization fraction and level populations for a certain electron temperature
	void calculateDensities(double T);

	double _n;
	const std::vector<double>& _isrf, _wavelength;

	double _ionizedFraction;

	TwoLevel _levels;
};

#endif /* _GASSPECIES_H_ */
