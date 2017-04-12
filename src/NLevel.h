#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"

#include <Eigen/Dense>
#include <array>
#include <map>
#include <string>

class NLevel
{
protected:
	/* All constants must be filled in by the derived class. The initializer list of this
	 * constructor will enable this by calling the functions below. The base class should
	 * implement these functions in such a way that the outputs are compatible. Note that these
	 * are not equivalent to getters, which are listed below them, as they generate the desired
	 * data from scratch. */
	NLevel(const Array& frequencyv, int nLv, const Eigen::VectorXd& ev,
	       const Eigen::VectorXd& gv, const Eigen::MatrixXd& avv,
	       const Eigen::MatrixXd& extraAvv);
	/* Another one, which does not directly set the frequency grid. setFrequencyv needs to be
	 * called after using this constructor. */
	NLevel(int nLv, const Eigen::VectorXd& ev, const Eigen::VectorXd& gv,
	       const Eigen::MatrixXd& avv, const Eigen::MatrixXd& extraAvv);
	virtual int makeNLv() const = 0;
	virtual Eigen::VectorXd makeEv() const = 0;
	virtual Eigen::VectorXd makeGv() const = 0;
	virtual Eigen::MatrixXd makeAvv() const = 0;
	virtual Eigen::MatrixXd makeExtraAvv() const = 0;

public:
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

public:
	/* Variables which are changed at every invocation of solveBalance or
	 * throughout the calculation, will be passed around using this struct
	 * (instead of storing them as members, which is not threadsafe. */
	typedef struct
	{
		double n, T;
		Eigen::MatrixXd bpvv, cvv;
		Eigen::VectorXd nv;
	} Solution;

	/* The matrices for this struct can be generated using the functions below */
private:
	/* Fill in the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the
	 Aij) and Pij is the line power (isrf integrated over the line profile) */
	Eigen::MatrixXd prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
	                                        const Eigen::MatrixXd& Cvv) const;

protected:
	/* Fill in the collision rates Cij. Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np (np
	 is the number of collision partners). Since these values depend on the model used, this
	 function has to be implemented in a derived class */
	virtual Eigen::MatrixXd prepareCollisionMatrix(double T, double electronDensity,
	                                               double protonDensity) const = 0;
	/* And the density vector nv corresponding to the level populations of the solution can be
	 * filled in using */
private:
	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji
	 + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term ne*np*alpha_i,
	 due to recombination. Returns the solution as a vector. */
	Eigen::VectorXd solveRateEquations(double n, const Eigen::MatrixXd& BPvv,
	                                   const Eigen::MatrixXd& Cvv,
	                                   const Eigen::VectorXd& sourceTerm,
	                                   const Eigen::VectorXd& sinkTerm,
	                                   int chooseConsvEq) const;

	/* Other functions */

	/* Abstraction of the loop over all lines. Executes thingWithLine for all combinations upper
	 * > lower that have _Avv(upper, lower) > 0. */
	void forAllLinesDo(std::function<void(size_t upper, size_t lower)> thingWithLine) const;

	double lineIntensityFactor(size_t upper, size_t lower, const Solution& info) const;
	double lineOpacityFactor(size_t upper, size_t lower, const Solution& info) const;

	/* Calculates the Voigt profile for a certain line, using the wavelengthgrid supplied at
	construction and the temperature and collision rates contained in the info struct. */
	Array lineProfile(size_t upper, size_t lower, const Solution& info) const;

	/* Or when the full solution is not yet known, and hence an Info object is not yet available
	 */
	Array lineProfile(size_t upper, size_t lower, double T, const Eigen::MatrixXd& Cvv) const;

	/* Calculates the contribution of A_ul to the total decay rate of u. This will determine the
	 probability that a photon is re-emitted. */
	double lineDecayFraction(size_t upper, size_t lower, const Solution& info) const;

	/* Variables which are the same for all invocations of solveBalance are stored as members */

	/* Wavelength grid */
	Array _frequencyv;
	/* Energy levels (constant) */
	int _nLv;
	Eigen::VectorXd _ev;
	/* Level degeneracy (constant) */
	Eigen::VectorXd _gv;
	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::MatrixXd _avv;
	/* Spontaneous transitions that do not produce line photons, but do influence the levels. A
	 * prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for single photon,
	 * but is about 8 s-1 for two photons*/
	Eigen::MatrixXd _extraAvv;

	/* A map allowing the naming of certain transitions. For every name, there can be multiple
	 * transitions contributing. For every transition there is one pair of level indices. */
	// std::map<std::string, std::vector<std::array<size_t, 2>>> _namedTransitions;

public:
	Array frequencyv() const { return _frequencyv; }
	void setFrequencyv(const Array& frequencyv) { _frequencyv = frequencyv; }

	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	Solution solveBalance(double atomDensity, double electronDensity, double protonDensity,
	                      double T, const Array& specificIntensityv, const Array& sourcev,
	                      const Array& sinkv) const;

	Array emissivityv(const Solution& info) const;
	Array opacityv(const Solution& info) const;
	Array scatteringOpacityv(const Solution& info) const;
};

#endif /* _NLEVEL_H_ */
