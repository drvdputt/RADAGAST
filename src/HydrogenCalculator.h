#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "GasState.h"
#include "Testing.h"
#include <vector>

class TwoLevel;
class FreeBound;
class FreeFree;

class HydrogenCalculator
{
public:
	/* Creates an object which can calculate the NLTE state, of a pure (ionized and atomic) hydrogen gas
	 and the resulting opacity and emission on the provided frequency grid. */
	HydrogenCalculator(const std::vector<double>& frequencyv);

	~HydrogenCalculator();

	/* Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature guess, and
	 a vector containing the radiation field in specific intensity per frequency units (on the same
	 frequency grid as the one provided at construction). */
	void solveBalance(double n, double Tinit, const std::vector<double>& specificIntensity);

	/* Exports the state of the gas as a compact, opaque object. Codes which make use of the gas module
	can use these objects to store a number of gas states. They can repeatedly give objects of this type
	to a single instance of the HydrogenCalculator to calculate any optical properties on-the-fly. The
	exact contents of a GasState and the way the optical properties are calculated (derived from densities
	vs caching them for example) are entirely up to the implementations of the functions below and the
	definition in GasState.h. */
	GasState exportState() const;

	/* The functions below hould provide a fast implementation to obtain the optical properties. The
	implementation will depend on what is stored in a GasState object. A good balance between the size of
	the GasState objects and the computation time needed for the optical properties needs to be found. */

	// 1 erg / cm3 = 0.1 J / m3
	double emissivity_SI(const GasState& gs, size_t iFreq)
	{
		return 0.1 * gs._emissivityv[iFreq];
	}
	// 1 / cm = 100 / m
	double opacity_SI(const GasState& gs, size_t iFreq)
	{
		return 100 * gs._opacityv[iFreq];
	}
	double scatteringOpacity_SI(const GasState& gs, size_t iFreq)
	{
		return 100 * gs._scatteringOpacityv[iFreq];
	}
	double absorptionOpacity_SI(const GasState& gs, size_t iFreq)
	{
		return 100 * (gs._opacityv[iFreq] - gs._scatteringOpacityv[iFreq]);
	}
	/* Finds a new balance, using information stored in the GasState to speed up the calculation. */
	void solveBalance(const GasState&, double n, double Tinit,
			const std::vector<double>& specificIntensity);

	/* Used by the balance solver to calculate the ionization fraction and level populations for a certain
	electron temperature, under influence of a blackbody isrf of that same temperature. Can be used by the client
	to manually set the temperature and calculate some properties which can be used as an initial guess. */
	void solveInitialGuess(double n, double T);


private:
	/* Give this testing function access to private members. */
	friend void Testing::testHydrogenCalculator();

	/* The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	std::vector<double> emissivityv() const;
	/* The total opacity at each frequency in 1 / cm */
	std::vector<double> opacityv() const;
	/* The scattering opacity used to simulate re-emission of line photons,
	 such as resonant scattering. The absorption opacity equals the total opacity minus this value. */
	std::vector<double> scatteringOpacityv() const;
	/* The amount of radiation scattered away by the scattering approach (in the same units as the
	emissivity) */
	std::vector<double> scatteredv() const;

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

	/* The bolometric emission by the continuum only (= cooling by recombination continuum) */
	double continuumEmission() const;
	/* The bolometric absorption by the continuum only (= ionization heating) */
	double continuumAbsorption() const;

	double np_ne() const
	{
		double ne = _ionizedFraction * _n;
		return ne*ne;
	}

	double nAtomic() const
	{
		return _n * (1. - _ionizedFraction);
	}

	void testHeatingCurve();

	void calculateDensities(double T);

private:
	/* To be set in constructor */
	const std::vector<double>& _frequencyv;

	/* To be set on invocation of solveBalance() */
	double _n{0};
	const std::vector<double>* _p_specificIntensityv{nullptr};

	/* Results of solveBalance() */
	double _T;
	double _ionizedFraction;

	/* Hide the other parts of the implementation from code that includes this header
	 to avoid the chaining of dependencies. */
	std::unique_ptr<TwoLevel> _levels{nullptr};

	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound{nullptr};
	std::unique_ptr<FreeFree> _freeFree{nullptr};
};

#endif /* _HYDROGENCALCULATOR_H_ */
