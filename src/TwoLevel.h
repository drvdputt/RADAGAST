#ifndef _TWOLEVEL_H_
#define _TWOLEVEL_H_
#include <Eigen/Dense>
#include <vector>

class TwoLevel {
public:
	// Creates an object that represents a two-level (component of a) medium. The level population equilibrium
	// will be calculated using the wavelength grid supplied as arguments
	// of the constructor.
	TwoLevel(const std::vector<double>& wavelengthv);

	// Calculates the level populations for a certain electron temperature and isrf. When ionization comes into play,
	// the recombination and ionization rates will act as source and sink terms for the levels. Therefore,
	// external sources or sink rates can be passed using vectors containing one number for each level.
	void doLevels(double n, double nc, double T,
			const std::vector<double>& isrf, const std::vector<double>& source,
			const std::vector<double>& sink);

	// Useful for the thermal balance. Is much faster than the calculation of the full emissio spectrum.
	double bolometricEmission(size_t upper, size_t lower) const;

	// The values needed for a radiative transfer cycle
	// The emission coefficient j_nu (erg/cm3/s/Hz)
	std::vector<double> calculateEmission() const;
	// The opacity alpha_nu, equivalent to kappaRho for dust (cm-1)
	std::vector<double> calculateOpacity() const;

private:
	// Calculates the Voigt profile for a certain line, using the wavelengthgrid supplied at construction
	// and the current temperature and collision rates.
	Eigen::ArrayXd lineProfile(size_t upper, size_t lower) const;

	// Fill in the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from the Aij)
	// and Pij is the line power (isrf integrated over the line profile)
	void prepareAbsorptionMatrix(const std::vector<double>& isrf);

	// Fill in the collision rates Cij
	// Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np
	// (np is the number of collision partners
	void prepareCollisionMatrix();

	// Setup and return the rate matrix M_ij = A_ji + B_ji * P_ji + C_ji.
	// Set up F and b using M_ij and the external source term ne*np*alpha_i, due to recombination
	// Stores the solution in the vector containing the densities _ni;
	void solveRateEquations(Eigen::Vector2d sourceTerm,
			Eigen::Vector2d sinkTerm, size_t chooseConsvEq);

	///////////////
	// Data members
	///////////////

	//--------------------------------
	// Constants, set at construction:
	//--------------------------------

	// Wavelength grid
	const std::vector<double>& _wavelengthv;
	// Energy levels (constant)
	Eigen::Vector2d _Ev;
	// Level degeneracy (constant)
	Eigen::Vector2d _gv;
	// A matrix (constant, lower triangle, zero diagonal)
	Eigen::Matrix2d _Avv;

	//----------------------------------------------------------
	// To be set/calculated at every iteration
	//----------------------------------------------------------
	// (i.e. everything that can change during the simulation / depends on the cell)
	// Total density
	double _n;
	// Density of collision partners (in this simple model, all partners are treated the same)
	double _nc;
	// Electron temperature
	double _T;
	// Line-averaged intensity for each pair of levels (line shape is temperature dependent!)
	Eigen::Matrix2d _BPvv;
	// Collisional transition rates
	Eigen::Matrix2d _Cvv;

	//-------------------------------------------------
	// Level densities (to be calculated in doLevels())
	//-------------------------------------------------
	Eigen::Vector2d _nv;
};

#endif
