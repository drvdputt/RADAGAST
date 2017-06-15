#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"
#include "LevelDataProvider.h"

#include <Eigen/Dense>
#include <array>
#include <map>
#include <memory>
#include <string>

class NLevel
{
protected:
	NLevel(LevelDataProvider* ldp);

public:
	virtual ~NLevel();

	int nLv() const { return _nLv; }

protected:
	Eigen::VectorXd ev() const { return _ev; }
	double ev(size_t i) const { return _ev(i); }

	Eigen::VectorXd gv() const { return _gv; }
	double gv(size_t i) const { return _gv(i); }

	Eigen::MatrixXd avv() const { return _avv; }
	double avv(size_t upper, size_t lower) { return _avv(upper, lower); }

	Eigen::MatrixXd extraAvv() const { return _extraAvv; }
	double extraAvv(size_t upper, size_t lower) const { return _extraAvv(upper, lower); }

	LevelDataProvider* ldp() const;

public:
	Array frequencyv() const { return _frequencyv; }
	void setFrequencyv(const Array& frequencyv) { _frequencyv = frequencyv; }

	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	/* Variables which are changed at every invocation of solveBalance or
	 throughout the calculation, will be passed around using this struct
	 (instead of storing them as members, which is not threadsafe. */
	typedef struct
	{
		/* The density and temperature for which the solution was calculated */
		double n, T;

		/* The density of each level population (cm-3) */
		Eigen::VectorXd nv;

		/* The induced radiative transition rates and collisional transition rates for this
		 configuratiom */
		Eigen::MatrixXd bpvv, cvv;
	} Solution;

	/* Calculates the level populations for a certain electron temperature and isrf. The number
	 of electrons and protons is needed to determine the collisional transition rates. When
	 ionization comes into play, the recombination and ionization rates will act as source and
	 sink terms for the levels. Therefore, external sources or sink rates can be passed using
	 vectors containing one number for each level. A struct of the type Solution as defined
	 above is returned */
	Solution solveBalance(double atomDensity, double electronDensity, double protonDensity,
	                      double T, const Array& specificIntensityv, const Array& sourcev,
	                      const Array& sinkv) const;

	/* The spectrum emitted by the line transitions, expressed as the emission coefficient j_nu
	f * (erg/cm3/s/Hz) */
	Array emissivityv(const Solution& s) const;

	/* The contribution to the emissivity non-line bound-bound processes, for example the
	 two-photon continuum 2s->1s transition of HI */
	virtual Array boundBoundContinuum(const Solution& s) const;

	/* The opacity alpha_nu, equivalent to kappaRho for dust (cm-1) */
	Array opacityv(const Solution& s) const;

	/* Heating rate due to collisional de-excitation (ergs / s / cm3) */
	double heating(const Solution& s) const;

	/* Cooling rate due to collisional excitation */
	double cooling(const Solution& s) const;

private:
	/* Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the
	 Aij) and Pij is the line power (isrf integrated over the line profile) */
	Eigen::MatrixXd prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
	                                        const Eigen::MatrixXd& Cvv) const;

protected:
private:
	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji
	 + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term ne*np*alpha_i,
	 due to recombination. Returns the solution as a vector. */
	Eigen::VectorXd solveRateEquations(double n, const Eigen::MatrixXd& BPvv,
	                                   const Eigen::MatrixXd& Cvv,
	                                   const Eigen::VectorXd& sourceTerm,
	                                   const Eigen::VectorXd& sinkTerm,
	                                   int chooseConsvEq) const;

	/* Abstraction of the loop over all lines. Executes thingWithLine for all combinations upper
	 > lower that have _Avv(upper, lower) > 0. If the levels are sorted, and all downward
	 transitions have line activity, then this function will loop over all elements of the
	 lower triangle of the level matrix. */
	void forActiveLinesDo(std::function<void(size_t upper, size_t lower)> thingWithLine) const;

	/* Calculates the intensity of a specific line [erg / s / cm2]. Multiplying with the line
	 profile will yield the specific intensity. */
	double lineIntensityFactor(size_t upper, size_t lower, const Solution& s) const;

	/* Computes the opacity of a line, not yet multiplied with the line profile */
	double lineOpacityFactor(size_t upper, size_t lower, const Solution& s) const;

	/* Calculates the Voigt profile for a certain line, using the wavelengthgrid supplied at
	construction and the temperature and collision rates contained in the Solution struct. */
	Array lineProfile(size_t upper, size_t lower, const Solution& s) const;

	/* Or when the full solution is not yet known, and hence a Solution object is not yet
	 available. */
	Array lineProfile(size_t upper, size_t lower, double T, const Eigen::MatrixXd& Cvv) const;

	/* Variables which are the same for all invocations of solveBalance are stored as members */

	/* Wavelength grid */
	Array _frequencyv;

	/* Energy levels (constant) */
	int _nLv{0};
	Eigen::VectorXd _ev;

	/* Level degeneracy (constant) */
	Eigen::VectorXd _gv;

	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::MatrixXd _avv;

	/* Spontaneous transitions that do not produce line photons, but do influence the levels. A
	 prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for single photon,
	 but is about 8 s-1 for two photons*/
	Eigen::MatrixXd _extraAvv;

	/* A polymorphic LevelDataProvider. The specific subclass that this data member is
	 initialized with depends on the subclass */
	std::unique_ptr<LevelDataProvider> _ldp;
};

#endif /* _NLEVEL_H_ */
