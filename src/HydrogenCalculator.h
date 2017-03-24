#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.h"
#include "GasState.h"

class TwoLevel;
class FreeBound;
class FreeFree;

class HydrogenCalculator
{
public:
	/* Creates an object which can calculate the NLTE state, of a pure (ionized and atomic) hydrogen gas
	 and the resulting opacity and emission on the provided frequency grid. */
	HydrogenCalculator(const Array& frequencyv);

	~HydrogenCalculator();

	/* Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature guess, and
	 a vector containing the radiation field in specific intensity per frequency units (on the same
	 frequency grid as the one provided at construction). */
	void solveBalance(double n, double Tinit, const Array& specificIntensity);

	/* Finds a new balance, using information stored in the GasState to speed up the calculation. */
	void solveBalance(const GasState&, double n, double Tinit, const Array& specificIntensity);

	/* Used by the balance solver to calculate the ionization fraction and level populations for a certain
	 electron temperature, under influence of a blackbody isrf of that same temperature. Can be used by the client
	 to manually set the temperature and calculate some properties which can be used as an initial guess. */
	void solveInitialGuess(double n, double T);

private:
	/* Calculates all the densities for a fixed temperature. Is repeatedly called by this class. */
	void calculateDensities(double T);

public:

	/* The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emissivityv() const;
	/* The total opacity at each frequency in 1 / cm */
	Array opacityv() const;
	/* The scattering opacity used to simulate re-emission of line photons,
	 such as resonant scattering. The absorption opacity equals the total opacity minus this value. */
	Array scatteringOpacityv() const;
	/* The amount of radiation scattered away by the scattering approach (in the same units as the
	 emissivity) */
	Array scatteredv() const;

	/* The total bolometric emission, in erg / s / cm^3, obtained by integrating the emissivity. */
	double emission() const;
	/* The total bolometric absorption, in erg / s / cm^3. This is an integral of the opacity
	 times the radiation field. */
	double absorption() const;

	/* The bolometric emission by the lines only. The emissivity of the photon re-emissions is also
	 included in this value. */
	double lineEmission() const;
	/* The bolometric absorption by the lines only. The absorption of re-emitted line photons is also
	 included here. */
	double lineAbsorption() const;
	/*  Taking the difference of the above two terms will cancel out the contributions of the
	 "scattered" photons and yield the heating/cooling contribution by the lines. */

	/* The bolometric emission by the continuum only (= cooling by recombination + free-free continuum) */
	double continuumEmission() const;
	/* The bolometric absorption by the continuum only (= ionization + free-free heating) */
	double continuumAbsorption() const;

	double temperature()
	{
		return _temperature;
	}

	double ionizedFraction()
	{
		return _ionizedFraction;
	}

	double np_ne() const
	{
		double ne = _ionizedFraction * _n;
		return ne * ne;
	}

	double nAtomic() const
	{
		return _n * (1. - _ionizedFraction);
	}

	void testHeatingCurve();

private:
	/* To be set in constructor */
	const Array& _frequencyv;

	/* To be set on invocation of solveBalance() */
	double _n
	{ 0 };
	const Array* _p_specificIntensityv
	{ nullptr };

	/* Results of solveBalance() */
	double _temperature;
	double _ionizedFraction;

	/* Hide the other parts of the implementation from code that includes this header
	 to avoid the chaining of dependencies. */
	std::unique_ptr<TwoLevel> _levels
	{ nullptr };

	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound
	{ nullptr };
	std::unique_ptr<FreeFree> _freeFree
	{ nullptr };
};

#endif /* _HYDROGENCALCULATOR_H_ */
