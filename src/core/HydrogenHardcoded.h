#ifndef GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_
#define GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_

#include "HydrogenDataProvider.h"

/** This provides data for a hydrogen level model in a hardcoded way. It is mainly for educational
    and compararive purposes. I used it to check if the implementation of \c HydrogenFromFiles
    correctly processed the j-resolved coefficients it read in, for example. */
class HydrogenHardcoded : public HydrogenDataProvider
{
public:
	HydrogenHardcoded();

	/** Returns the number of levels */
	size_t numLv() const override;

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

	int nMax() const override { return 5; }

	size_t indexOutput(int n, int l) const override;

	std::array<size_t, 2> twoPhotonIndices() const override;

	//-----------------------------------//
	// FUNCTIONS RETURNING VARIABLE DATA //
	//-----------------------------------//
	/* These functions provide coefficients that depend on external variables such as the
	   temperature. */

	/** Returns a matrix containing the collisional transition rates (already multiplied with
	    the partner density), for a given temperature and proton and electron densities. */
	EMatrix cvv(const GasStruct& gas) const override;

private:
	EVector the_ev;
	EVector the_gv;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENHARDCODED_H_ */
