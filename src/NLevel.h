#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"
#include "LevelDataProvider.h"

#include <Eigen/Dense>
#include <array>
#include <map>
#include <memory>
#include <string>

/* This class contains a generic implementation for calculating the statistical equilibrium of a
   system with a certain amount of levels (private part), and some properties that result from it
   (public part).

   It contains data members to store all the constants needed for these calculations: The number of
   levels, the energies of these levels and their degeneracies, and the spontaneous transition rates
   between them. Since an object of this clas is to be reused many times, or by multiple threads at
   the same time, this class is designed to be stateless: notice how almost all of the methods are
   const.

   However, the calculation of the statisitical equilibrium does need some non-constant variables
   such as the temperature, population densities or the current collision rates. These are stored
   into an object of the type NLevel::Solution, a struct defined in this class. A Solution object
   can be obtained by calling the solveBalance function. By giving this object to one of the public
   functions (e.g. emissivityv), various properties that depend on the densities, temperature,
   collision rates... can be calculated.

   All of the necessary data for the equilibrium calculation is obtained through a class called
   LevelDataProvider. The latter is an abstract class, and the type of atom/molecule that is
   simulated by the implementation in NLevel class will depend on the subclass of LevelDataProvider
   that is used. Each subclass of NLevel will determine what kind of LevelDataProvider is used
   during construction. Some functions of LevelDataProvider are called during construction of the
   NLevel object to set the various constants. In contrast, the functions of the LevelDataProvider
   for the collision coefficients are called at every invocation of solveBalance, as these depend on
   the temperature and collision partner densities.

   Subclasses are allowed to override some of the implementations in this class, but they must use
   these same data members. That's why a set of protected getters is provided. */
class NLevel
{
protected:
	/* A subclass will create a new LevelDataProvider object, and provide a pointer to it to the
	   NLevel constructor. The number of levels, energies, multiplicities, and transition
	   coefficients are filled in. */
	NLevel(LevelDataProvider* ldp);

public:
	virtual ~NLevel();

	/* Return the number of levels in the solution */
	int numLv() const { return _numLv; }

protected:
	/* A number of protected getters are provided, so the subclasses can make use of these
	   coefficients. */
	Eigen::VectorXd ev() const { return _ev; }
	double ev(size_t i) const { return _ev(i); }

	Eigen::VectorXd gv() const { return _gv; }
	double gv(size_t i) const { return _gv(i); }

	Eigen::MatrixXd avv() const { return _avv; }
	double avv(size_t upper, size_t lower) { return _avv(upper, lower); }

	Eigen::MatrixXd extraAvv() const { return _extraAvv; }
	double extraAvv(size_t upper, size_t lower) const { return _extraAvv(upper, lower); }

public:
	/* Getter and setter for the frequency grid to perform the calculations on. All of the input
	   (ouput) spectra must have (will have) the same frequency points. */
	Array frequencyv() const { return _frequencyv; }
	void setFrequencyv(const Array& frequencyv) { _frequencyv = frequencyv; }

	/* Ouputs some properties about the different line transitions taken into account by this
	   instance of NLevel. The results for the number of lines, their frequencies and their
	   natural width are returned by reference. */
	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	/* Variables which are changed at every invocation of solveBalance or throughout the
	   calculation, will be passed around using this struct (instead of storing them as members,
	   which is not threadsafe. Note that Solution objects are not interchangeable between
	   NLevel instances, as the number of levels and their properties can be different.

	   TODO: Reconsider how we can prevent breaking abstraction here. Ideally, we'd want a
	   mechanism so that Solutions can only be used with the object that created them. */
	typedef struct
	{
		/* The total density and temperature of the ensemble of atoms/molecules for which
		   the solution was calculated */
		double n, T;

		/* The density of each level population (cm-3) */
		Eigen::VectorXd nv;

		/* The induced radiative transition rates and collisional transition rates for this
		   configuration. These are needed to calculate for example the line broadening. */
		Eigen::MatrixXd bpvv, cvv;
	} Solution;

	/* Calculates the level populations for a certain electron temperature and isrf. The
	   temperature, and the electron and proton density are needed to determine the collisional
	   transition rates and recombination rates (this format will probably change when we want
	   to use the same class for other atoms or by extension molecules). A struct of the type
	   Solution is returned, with all of its members correctly filled in. The collision and
	   recombination data to perform the calculation are obtained from the LevelDataProvider,
	   while the induced transitions rates (B coefficients * line power) are derived from the
	   given specific intensity. */
	Solution solveBalance(double atomDensity, double electronDensity, double protonDensity,
	                      double T, const Array& specificIntensityv) const;

	/* The spectrum emitted by the line transitions, expressed as the emission coefficient j_nu
	   f * (erg/cm3/s/Hz). */
	Array emissivityv(const Solution& s) const;

	/* The contribution to the emissivity non-line bound-bound processes, for example the
	   two-photon continuum 2s->1s transition of HI. */
	virtual Array boundBoundContinuum(const Solution& s) const;

	/* The opacity alpha_nu, equivalent to kappaRho for dust (cm-1). */
	Array opacityv(const Solution& s) const;

	/* Heating rate due to collisional de-excitation (ergs / s / cm3). */
	double heating(const Solution& s) const;

	/* Cooling rate due to collisional excitation. */
	double cooling(const Solution& s) const;

private:
	/* Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the
	   Aij) and Pij is the line power (isrf integrated over the line profile). */
	Eigen::MatrixXd prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
	                                        const Eigen::MatrixXd& Cvv) const;

private:
	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji
	   + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term
	   ne*np*alpha_i, due to recombination. Returns the solution as a vector. */
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
	int _numLv{0};
	Eigen::VectorXd _ev;

	/* Level degeneracy (constant) */
	Eigen::VectorXd _gv;

	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::MatrixXd _avv;

	/* Spontaneous transitions that do not produce line photons, but do influence the levels. A
	   prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for single photon,
	   but is about 8 s-1 for two photons. */
	Eigen::MatrixXd _extraAvv;

	/* A polymorphic LevelDataProvider. The specific subclass that this data member is
	   initialized with depends on the subclass. */
	std::unique_ptr<LevelDataProvider> _ldp;
};

#endif /* _NLEVEL_H_ */
