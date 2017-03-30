#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"

#include <Eigen/Dense>

class NLevel
{
public:
	NLevel(const Array& frequencyv);

	int N() {return _nLv;}

	void solveBalance(double atomDensity, double electronDensity, double protonDensity,
	                  double T, const Array& specificIntensityv, const Array& sourcev,
	                  const Array& sinkv);

	Array emissivityv() const;

	Array opacityv() const;

	Array scatteringOpacityv() const;
	Array absorptionOpacityv() const;

private:
	double lineIntensityFactor(size_t upper, size_t lower) const;
	double lineOpacityFactor(size_t upper, size_t lower) const;

	/* Calculates the Voigt profile for a certain line, using the wavelengthgrid supplied at
	construction and the current temperature and collision rates. */
	Array lineProfile(size_t upper, size_t lower) const;

	/* Calculates the contribution of A_ul to the total decay rate of u. This will determine the
	 probability that a photon is re-emitted. */
	double lineDecayFraction(size_t upper, size_t lower) const;

	/* Fill in the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the
	 Aij) and Pij is the line power (isrf integrated over the line profile) */
	void prepareAbsorptionMatrix(const Array& specificIntensityv);

	void setElectronCollisionRate(size_t upper, size_t lower, double bigGamma);

	/* Fill in the collision rates Cij. Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np (np
	 is the number of collision partners). */
	void prepareCollisionMatrix();

	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji
	 + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term ne*np*alpha_i,
	 due to recombination. Stores the solution in the vector containing the densities _ni. */
	void solveRateEquations(Eigen::VectorXd sourceTerm, Eigen::VectorXd sinkTerm,
	                        int chooseConsvEq);

	/* Wavelength grid */
	const Array& _frequencyv;
	/* Energy levels (constant) */
	int _nLv{6};
	Eigen::VectorXd _Ev;
	/* Level degeneracy (constant) */
	Eigen::VectorXd _gv;
	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::MatrixXd _Avv;
	/* Spontaneous transitions that do not produce line photons, but do influence the levels. A
	 * prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for single photon,
	 * but is about 8 s-1 for two photons*/
	Eigen::MatrixXd _extraAvv;

	/* To be set or calculated at every invocation of the balance calculation (i.e. everything
	 that can change during the simulation or depends on the cell) */

	/* Total density */
	double _n{0};
	/* Density of collision partners, in this case electrons and protons */
	double _ne{0};
	double _np{0};
	/* Electron temperature */
	double _T{0};
	/* Line-averaged intensity for each pair of levels (line shape is temperature dependent, so
	 this has to be updated whenever it changes) */
	Eigen::MatrixXd _BPvv;
	/* Collisional transition rates */
	Eigen::MatrixXd _Cvv;

	/* Level densities (to be calculated) */
	Eigen::VectorXd _nv;
};

#endif /* _NLEVEL_H_ */
