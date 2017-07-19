#ifndef GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_
#define GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_

#include "HydrogenDataProvider.h"

class HydrogenHardcoded : public HydrogenDataProvider
{
public:
	HydrogenHardcoded();

	/** Returns the number of levels */
	int numLv() const override;

	/** Returns a vector containing the energy of each level */
	EVector ev() const override;
	/** Returns a vector containing the degeneracy of each level */
	EVector gv() const override;

	/** Returns a matrix containing the Einstein A coefficients for all levels. Indexed on
	    (upper, lower), making it a lower triangle matrix. */
	EMatrix avv() const override;

	/** Returns a matrix containing any extra spontaneous decays between levels. This matrix can
	    be used to describe spontaneous decays that do NOT produce line radiation (for example
	    two-photon processes, which generate a continuum instead). */
	EMatrix extraAvv() const override;

	std::array<int, 2> twoPhotonIndices() const override;

	//-----------------------------------//
	// FUNCTIONS RETURNING VARIABLE DATA //
	//-----------------------------------//
	/* These functions provide coefficients that depend on external variables such as the
	   temperature. */

	/** Returns a matrix containing the collisional transition rates (already multiplied with
	    the partner density), for a given temperature and proton and electron densities. */
	EMatrix cvv(double T, double ne, double np) const override;

	/** Returns a vector containing partial recombination rates based on fits I found somewhere
	    (see source code for origin). */
	EVector sourcev(double T, double ne, double np) const override;
	EVector sinkv(double T, double ne, double np) const override;

private:
	EVector the_ev;
	EVector the_gv;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_ */
