#ifndef GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_
#define GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_

#include "LevelDataProvider.h"

class HydrogenHardcoded : public LevelDataProvider
{
public:
	HydrogenHardcoded();

	/* Returns the number of levels */
	int numLv() const override;
	/* Returns a vector containing the energy of each level */
	Eigen::VectorXd ev() const override;
	/* Returns a vector containing the degeneracy of each level */
	Eigen::VectorXd gv() const override;

	/* Returns a matrix containing the Einstein A coefficients for all levels. Indexed on
	   (upper, lower), making it a lower triangle matrix. */
	Eigen::MatrixXd avv() const override;

	/* Returns a matrix containing any extra spontaneous decays between levels. This matrix can
	   be used to describe spontaneous decays that do NOT produce line radiation (for example
	   two-photon processes, which generate a continuum instead). */
	Eigen::MatrixXd extraAvv() const override;

	//-----------------------------------//
	// FUNCTIONS RETURNING VARIABLE DATA //
	//-----------------------------------//
	/* These functions provide coefficients that depend on external variables such as the
	   temperature. */

	/* Returns a matrix containing the collisional transition rates (already multiplied with the
	   partner density), for a given temperature and proton and electron densities. */
	Eigen::MatrixXd cvv(double T, double ne, double np) const override;

	/* Returns a vector containing partial recombination rates based on fits I found somewhere
	   (see source code for origin). */
	Eigen::VectorXd alphav(double T) const override;

private:
	Eigen::VectorXd the_ev;
	Eigen::VectorXd the_gv;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_ */
