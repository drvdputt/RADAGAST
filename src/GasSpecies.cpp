#include "GasSpecies.h"
#include "IonizationBalance.h"
#include <iostream>
#include <vector>

GasSpecies::GasSpecies(double n, const std::vector<double>& isrf, const std::vector<double>& wavelength)
	: _n(n), _isrf(isrf), _wavelength(wavelength), _ionizedFraction(0), _levels(wavelength, isrf)
{
}

void GasSpecies::solveBalance()
{
	// Initial guess for the temperature
	double T = 1000;

	calculateDensities(T);
}

std::vector<double> GasSpecies::emissivity() const
{


	return _levels.calculateEmission();
}

std::vector<double> GasSpecies::opacity() const
{
	return _levels.calculateOpacity();
}

void GasSpecies::calculateDensities(double T)
{
	_ionizedFraction = Ionization::ionizedFraction(_n, T, _wavelength, _isrf);

	std::cout << "Ionized fraction = " << _ionizedFraction << std::endl;

	double nAtomic = _n * (1 - _ionizedFraction);
	// For now we can sum over the collision partners, as TwoLevel threats all collision parnters in the same way
	// neutral + ion + electron densities
	double nTotal = (1 + _ionizedFraction) * _n;
	double np_ne_alpha = _n*_n*_ionizedFraction*_ionizedFraction * Ionization::recombinationRate(T);
	_levels.doLevels(nAtomic, nTotal, T, np_ne_alpha);
}
