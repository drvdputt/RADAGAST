#ifndef _GASSPECIES_H_
#define _GASSPECIES_H_

#include "TwoLevel.h"

#include <vector>

class HydrogenCalculator
{
public:
	// Creates an object which can calculate the NLTE state, of a pure (ionized and atomic) hydrogen gas
	// and the resulting opacity and emission on the provided frequency grid.
	HydrogenCalculator(const std::vector<double>& frequencyv);

	// Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature guess, and
	// a vector containing the radiation field in specific intensity per frequency units (on the same frequency grid
	// as the one provided at construction).
	void solveBalance(double n, double Tinit, const std::vector<double>& specificIntensity);

	// The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz.
	std::vector<double> emissivityv() const;
	// The total opacity at each frequency in 1 / cm
	std::vector<double> opacityv() const;
	// The scattering opacity used to simulate re-emission of line photons,
	// such as resonant scattering.
	// The absorption opacity equals the total opacity minus this value.
	std::vector<double> scatteringOpacityv() const;

	// The total bolometric emission, in erg / s / cm^3, obtained by integrating the emissivity.
	double emission() const;
	// The total bolometric absorption, in erg / s / cm^3. This is an integral of the opacity
	// times the radiation field.
	double absorption() const;

	// The bolometric emission by the lines only. The emissivity of the photon re-emissions is also included in this value.
	double lineEmission() const;
	// The bolometric absorption by the lines only. The absorption of re-emitted line photons is also
	// included here.
	double lineAbsorption() const;
	// Taking the difference of these two terms will cancel out the contributions of the
	// "scattered" photons and yield the heating/cooling contribution by the lines.

	// The bolometric emission by the continuum only (= cooling by recombination continuum)
	double continuumEmission() const;
	// The bolometric absorption by the continuum only (= ionization heating)
	double continuumAbsorption() const;

	void testHeatingCurve();

private:
	// Calculates the ionization fraction and level populations for a certain electron temperature and isrf
	void calculateDensities(double T);

	// Calculate the emission coefficient for the optical recombination continuum, interpolated from data.
	// Returned in units [density^-1][power]/[frequency interval] cm^3 erg / s / cm.
	// The emissivity ([power][density]/[frequency interval]) can be obtained by multiplying this value with ne_np
	std::vector<double> recombinationEmissionCoeffv(double T) const;

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
