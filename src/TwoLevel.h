#ifndef _TWOLEVEL_H_
#define _TWOLEVEL_H_

#include "Array.h"

#include <Eigen/Dense>

class TwoLevel
{
public:
	/* Creates an object that represents a two-level (component of a) medium. The level population
	 equilibrium will be calculated using the frequency grid supplied as argument of the constructor. */
	TwoLevel(const Array& frequencyv);

	/* Calculates the level populations for a certain electron temperature and isrf. The number of
	 electrons and protons is needed to determine the collisional transition rates. When ionization comes
	 into play, the recombination and ionization rates will act as source and sink terms for the levels.
	 Therefore, external sources or sink rates can be passed using vectors containing one number for each
	 level. */
	void solveBalance(double n, double ne, double np, double T, const Array& specificIntensity,
			const Array& source, const Array& sink);

	/* Useful for the thermal balance. Is much faster than the calculation of the full emission spectrum.
	 */
	double lineIntensity(size_t upper, size_t lower) const;

	/* The values needed for a radiative transfer cycle */

	/* The emission coefficient j_nu (erg/cm3/s/Hz) */
	Array totalEmissivityv() const;

	/* The opacity alpha_nu, equivalent to kappaRho for dust (cm-1) */
	Array totalOpacityv() const;

	/* The part of the opacity which acts as a source of scattering. This scattering is equivalent to the
	 immediate emission of a photon by the line that just absorbed it. */
	Array scatteringOpacityv() const;
	Array absorptionOpacityv() const;

private:
	/* Calculates the Voigt profile for a certain line, using the wavelengthgrid supplied at construction
	 and the current temperature and collision rates. */
	Eigen::ArrayXd lineProfile(size_t upper, size_t lower) const;

	/* Calculates the contribution of A_ul to the total decay rate. This will determine the probability
	 that a photon is re-emitted. */
	double radiativeDecayFraction(size_t upper, size_t lower) const;

	/* Fill in the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the Aij)
	 and Pij is the line power (isrf integrated over the line profile) */
	void prepareAbsorptionMatrix(const Array& specificIntensity);

	/* Fill in the collision rates Cij. Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np (np is the
	 number of collision partners). */
	void prepareCollisionMatrix();

	/* Following the notation of the gasPhysics document, construct the rate matrix M_ij = A_ji + B_ji *
	 P_ji + C_ji. Set up F and b using M_ij and the external source term ne*np*alpha_i, due to
	 recombination. Stores the solution in the vector containing the densities _ni. */
	void solveRateEquations(Eigen::Vector2d sourceTerm, Eigen::Vector2d sinkTerm, size_t chooseConsvEq);

	/* Wavelength grid */
	const Array& _frequencyv;
	/* Energy levels (constant) */
	Eigen::Vector2d _Ev;
	/* Level degeneracy (constant) */
	Eigen::Vector2d _gv;
	/* A matrix (constant, lower triangle, zero diagonal) */
	Eigen::Matrix2d _Avv;

	/* To be set or calculated at every invocation of the balance calculation (i.e. everything that can
	 change during the simulation or depends on the cell) */

	/* Total density */
	double _n;
	/* Density of collision partners, in this case electrons and protons */
	double _ne;
	double _np;
	/* Electron temperature */
	double _T;
	/* Line-averaged intensity for each pair of levels (line shape is temperature dependent, so this has
	 to be updated whenever it changes) */
	Eigen::Matrix2d _BPvv;
	/* Collisional transition rates */
	Eigen::Matrix2d _Cvv;

	/* Level densities (to be calculated) */
	Eigen::Vector2d _nv;
};

#endif /* _TWOLEVEL_H_ */
