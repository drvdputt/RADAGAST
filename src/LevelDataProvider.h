#ifndef _SRC_LEVELDATAPROVIDER_H_
#define _SRC_LEVELDATAPROVIDER_H_

#include <Eigen/Dense>

/* When an object of the class NLevel is set up, it needs some constant data, and some data that
can vary with temperature. This data often needs to be derived from other data that was read in from
files. To abstract the read-in of files, and the conversion of their data to a usable format for the
calculations, this class has been created. Subclasses of
LevelDataProvider should read in a set of data, or have some data hardcoded, and provide
implementations of the functions below. NLevel can then call these functions retrieve data that can
be used for the statistical equilibrium calculation.  Subclasses may provide additional
configuration steps
which determine the actual number of levels. For example, one could configure a certain number of
l-resolved levels, and have the rest collapsed. Hence, the data that
is returned by these functions can be of a varying dimension (numbe of levels) or content (the
levels/terms treated can be different)
for the same
subclass.
*/
class LevelDataProvider
{
protected:
	/* During construction, a subclass will read in a certain dataset, representing a system
	 with line transitions. */
	LevelDataProvider();

public:
	virtual ~LevelDataProvider();

	//-----------------------------------//
	// FUNCTIONS RETURNING CONSTANT DATA //
	//-----------------------------------//
	/* These functions convert the internally loaded data to usable matrices an vectors for
	 level calculations. These are usually called only once, during the setup of a typical run.
	 */

	/* Returns the number of levels */
	virtual int numLv() const = 0;

	/* Returns a vector containing the energy of each level */
	virtual Eigen::VectorXd ev() const = 0;

	/* Returns a vector containing the degeneracy of each level */
	virtual Eigen::VectorXd gv() const = 0;

	/* Returns a matrix containing the Einstein A coefficients for all levels. Indexed on
	 (upper,
	 lower), making it a lower triangle matrix. */
	virtual Eigen::MatrixXd avv() const = 0;

	/* Returns a matrix containing any extra spontaneous decays between levels. This matrix can
	 be used to describe spontaneous decays that do NOT produce line radiation (for example
	 two-photon processes, which generate a continuum instead). */
	virtual Eigen::MatrixXd extraAvv() const = 0;

	//-----------------------------------//
	// FUNCTIONS RETURNING VARIABLE DATA //
	//-----------------------------------//
	/* These functions provide coefficients that depend on external variables such as the
	 temperature */

	/* Returns a matrix containing the collisional transition rates (already multiplied with the
	 partner density), for a given temperature and proton and electron densities. Calculate the
	 collision rates Cij. Rij = Cij * ni = q_ij * np * ni --> Cij = q_ij * np
	 with np the density of collision partners */
	virtual Eigen::MatrixXd cvv(double T, double ne, double np) const = 0;
};

#endif /* _SRC_LEVELDATAPROVIDER_H_ */
