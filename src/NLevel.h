#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"

#include <Eigen/Dense>

class NLevel
{
public:
	/* Variables which are changed and used thoughout the calculation will be passed around
	 * using this struct (instead of setting them as members, which is not threadsafe. */
	typedef struct
	{
		double n, T;
		Eigen::MatrixXd BPvv, Cvv;
		Eigen::VectorXd nv;
	} Solution;

	NLevel();
	NLevel(const Array& frequencyv);

	int N() const { return _nLv; }
	Array frequencyv() const { return _frequencyv; }
	void setFrequencyv(const Array& frequencyv) { _frequencyv = frequencyv; }

	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	Solution solveBalance(double atomDensity, double electronDensity, double protonDensity,
	                  double T, const Array& specificIntensityv, const Array& sourcev,
	                  const Array& sinkv) const;

	Array emissivityv(const Solution& info) const;
	Array opacityv(const Solution& info) const;
	Array scatteringOpacityv(const Solution& info) const;

private:
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

	/* Fill in the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the
	 Aij) and Pij is the line power (isrf integrated over the line profile) */
	Eigen::MatrixXd prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
	                                        const Eigen::MatrixXd& Cvv) const;

	/* Fill in the collision rates Cij. Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np (np
	 is the number of collision partners). */
	Eigen::MatrixXd prepareCollisionMatrix(double T, double electronDensity,
	                                       double protonDensity) const;

	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji
	 + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term ne*np*alpha_i,
	 due to recombination. Returns the solution as a vector. */
	Eigen::VectorXd solveRateEquations(double n, const Eigen::MatrixXd& BPvv,
	                                   const Eigen::MatrixXd& Cvv,
	                                   const Eigen::VectorXd& sourceTerm,
	                                   const Eigen::VectorXd& sinkTerm,
	                                   int chooseConsvEq) const;

	/* Wavelength grid */
	Array _frequencyv;
	/* Energy levels (constant) */
	const int _nLv{6};
	Eigen::VectorXd _Ev;
	/* Level degeneracy (constant) */
	Eigen::VectorXd _gv;
	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::MatrixXd _Avv;
	/* Spontaneous transitions that do not produce line photons, but do influence the levels. A
	 * prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for single photon,
	 * but is about 8 s-1 for two photons*/
	Eigen::MatrixXd _extraAvv;
};

#endif /* _NLEVEL_H_ */
