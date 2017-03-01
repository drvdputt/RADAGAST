#ifndef _GASSPECIES_H_
#define _GASSPECIES_H_

#include "TwoLevel.h"

#include <vector>

class GasSpecies
{
public:
	GasSpecies(const std::vector<double>& frequencyv);

	void solveBalance(double n, double Tinit, const std::vector<double>& specificIntensity);

	std::vector<double> emissivityv() const;
	std::vector<double> opacityv() const;

	double emission() const;
	double absorption() const;

	double lineEmission() const;
	double lineAbsorption() const;

	double continuumEmission() const;
	double continuumAbsorption() const;

	void testHeatingCurve();

private:
	// Calculates the ionization fraction and level populations for a certain electron temperature and isrf
	void calculateDensities(double T);

	// Calculate the emission coefficient for the optical recombination continuum, interpolated from data.
	// Returned in units [density^-1][power]/[frequency interval] cm^3 erg / s / cm.
	// The emissivity ([power][density]/[frequency interval]) can be obtained by multiplying this value with ne_np
	std::vector<double> recombinationEmissionCoeff(double T) const;

	// To be set in constructor
	const std::vector<double>& _frequencyv;

	// Data to be loaded in constructor body
	// First index is for frequency, second for temperature
	std::vector<std::vector<double>> _gammaDaggervv;
	// Vector containing the threshold frequencies, i.e. those of the lines starting with 1.
	// Is needed for applying equation 1 of Ercolano and Storey 2006 (MNRAS 372, 1875)
	std::vector<double> _thresholdv;

	// The accompanying log-temperature grid
	std::vector<double> _logTemperaturev;

	// To be set on invocation of solveBalance()
	double _n;
	const std::vector<double>* _p_specificIntensityv;

	// Results of solveBalance()
	double _T;
	double _ionizedFraction;
	TwoLevel _levels;
};

#endif /* _GASSPECIES_H_ */
