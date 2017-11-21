#ifndef GASMODULE_GIT_SRC_LEVELDATAPROVIDER_H_
#define GASMODULE_GIT_SRC_LEVELDATAPROVIDER_H_

#include "EigenAliases.h"

/** When an object of the class NLevel is set up, it needs some constant data, and some data that
    can vary with temperature. This data often needs to be derived from other data that was read in
    from files. To abstract the read-in of files, and the conversion of their data to a usable
    format for the calculations, this class has been created. Subclasses of LevelDataProvider should
    read in a set of data, or have some data hardcoded, and provide implementations of the functions
    below. NLevel can then call these functions retrieve data that can be used for the statistical
    equilibrium calculation.  Subclasses may provide additional configuration steps which determine
    the actual number of levels. For example, one could configure a certain number of l-resolved
    levels, and have the rest collapsed. Hence, the data that is returned by these functions can be
    of a varying dimension (number of levels) or content (the levels/terms treated can be different)
    for the same subclass. */
class LevelDataProvider
{
protected:
	LevelDataProvider();

public:
	virtual ~LevelDataProvider();

	//-----------------------------------------------//
	// PUBLIC FUNCTIONS AKA THE OUTPUT OF THIS CLASS //
	//-----------------------------------------------//

	// FUNCTIONS RETURNING CONSTANT DATA //
	/* These functions convert the data stored in an object of this class to usable matrices an
	   vectors for level calculations. These are usually called only once, during the setup of a
	   typical run. */

	/** Returns the number of levels. */
	virtual size_t numLv() const = 0;

	/** Returns a vector containing the energy of each level. [erg]*/
	virtual EVector ev() const = 0;

	/** Returns a vector containing the degeneracy of each level. */
	virtual EVector gv() const = 0;

	/** Returns a matrix containing the Einstein A coefficients for all levels. Indexed on
	    (upper, lower), making it a lower triangle matrix. [s-1] */
	virtual EMatrix avv() const = 0;

	/** Returns a matrix containing any extra spontaneous decays between levels. This matrix can
	    be used to describe spontaneous decays that do NOT produce line radiation (for example
	    two-photon processes, which generate a continuum instead). [s-1] */
	virtual EMatrix extraAvv() const = 0;

	// FUNCTIONS RETURNING VARIABLE DATA //
	/** These functions provide coefficients that depend on external variables such as the
	    temperature. */

	/** Returns a matrix containing the collisional transition rates (already multiplied with
	    the partner density), for a given temperature and proton and electron densities.
	    Calculate the collision rates Cij. The rate we need, Rij = Cij * ni = q_ij * np * ni -->
	    Cij = q_ij * np with np the density of collision partner. This function return Cij in
	    these equations, hence the unit is [s-1]. */
	virtual EMatrix cvv(double T, double ne, double np) const = 0;

	/** Return a vector which describes the rate at which new atoms/molecules appear directly
	    into the levels, due to external processes such as recombination or the atom/molecule
	    being a reaction product. These values will appear as constants in the equilibrium
	    equations. [cm-3 s-1] */
	virtual EVector sourcev(double T, double ne, double np) const = 0;

	/** Returns a vector which describes the fractional rate at which atoms/molecules will
	    disappear from a level, for example because they are being ionized, or because they are
	    involed in a chemical reaction. These values typically depend on the radiation field and
	    other environmental parameters. In the equilibrium equations, this value will be
	    muliplied by the level population density, hence the unit needs to be [s-1]. */
	virtual EVector sinkv(double T, double n, double ne, double np) const = 0;
};

#endif /* GASMODULE_GIT_SRC_LEVELDATAPROVIDER_H_ */
