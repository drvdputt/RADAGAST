#ifndef _HYDROGENCALCULATOR_H_
#define _HYDROGENCALCULATOR_H_

#include "Array.h"
#include "GasState.h"
#include "NLevel.h"

#include <memory>

class FreeBound;
class FreeFree;

class GasInterfaceImpl
{
public:
	typedef struct
	{
		double n, T, f;
		Array specificIntensityv;
		NLevel::Solution levelSolution;
	} Solution;

	/* Creates an object which can calculate the NLTE state, of a pure (ionized and atomic)
	 hydrogen gas and the resulting opacity and emission on the provided frequency grid. */
	GasInterfaceImpl(const Array& frequencyv);

	/* Creates a HydrogenCalculator and creates a frequency grid / adjusts the given frequency
	 * grid (by reference) if suggestedGrid = true.*/
	GasInterfaceImpl(const Array& frequencyv, bool improveGrid);

	Array frequencyv() const { return _frequencyv; }

	~GasInterfaceImpl();

	/* Solves for the NLTE, given a total hydrogen density n, an initial (electron) temperature
	 guess, and a vector containing the radiation field in specific intensity per frequency
	 units (on the same frequency grid as the one provided at construction). */
	void solveBalance(GasState&, double n, double Tinit, const Array& specificIntensity) const;

	/* Used by the balance solver to calculate the ionization fraction and level populations for
	 a certain electron temperature, under influence of a blackbody isrf of that same
	 temperature. Can be used by the client to manually set the temperature and calculate some
	 properties which can be used as an initial guess. */
	void solveInitialGuess(GasState&, double n, double T) const;

private:
	/* Calculates all the densities for a fixed temperature. Is repeatedly called by this class.
	 */
	Solution calculateDensities(double n, double T, const Array& specificIntensityv) const;

public:
	/* The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emissivityv(const Solution&) const;
	/* The total opacity at each frequency in 1 / cm */
	Array opacityv(const Solution&) const;
	/* A possible scattering opacity */
	Array scatteringOpacityv(const Solution&) const;

	/* The total bolometric emission, in erg / s / cm^3, obtained by integrating the emissivity.
	 */
	double cooling(const Solution&) const;
	/* The total bolometric absorption, in erg / s / cm^3. This is an integral of the opacity
	 times the radiation field. */
	double heating(const Solution&) const;

	/* The bolometric emission by the lines only. The emissivity of the photon re-emissions is
	 also included in this value. */
	double lineCooling(const Solution&) const;
	/* The bolometric absorption by the lines only. The absorption of re-emitted line photons is
	 also included here. */
	double lineHeating(const Solution&) const;

	/* The bolometric emission by the continuum only (= cooling by recombination + free-free
	 * continuum) */
	double continuumCooling(const Solution&) const;
	/* The bolometric absorption by the continuum only (free-free heating) in erg / s / cm3 */
	double continuumHeating(const Solution&) const;

	inline double np_ne(const Solution& s) const
	{
		double np = s.n * s.f;
		return np * np;
	}

	inline double nAtomic(const Solution& s) const { return (1 - s.f) * s.n; }

	void testHeatingCurve(double n, const Array& specificIntensityv) const;

private:
	/* To be set in constructor */
	Array _frequencyv;

	/* Pointers to other parts of the implementation, to make late initialization possible */
	std::unique_ptr<NLevel> _boundBound;
	/* Continuum contributions */
	std::unique_ptr<FreeBound> _freeBound;
	std::unique_ptr<FreeFree> _freeFree;
};

#endif /* _HYDROGENCALCULATOR_H_ */
