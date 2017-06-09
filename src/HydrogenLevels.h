#ifndef _SRC_HYDROGENLEVELS_H_
#define _SRC_HYDROGENLEVELS_H_

#include "NLevel.h"

class HydrogenLevels : public NLevel
{
public:
	HydrogenLevels();
	HydrogenLevels(const Array& frequencyv);

protected:
	int makeNLv() const override;
	Eigen::VectorXd makeEv() const override;
	Eigen::VectorXd makeGv() const override;
	Eigen::MatrixXd makeAvv() const override;
	Eigen::MatrixXd makeExtraAvv() const override;

	Eigen::MatrixXd prepareCollisionMatrix(double T, double electronDensity,
	                                       double protonDensity) const override;
	Array boundBoundContinuum(const Solution& s) const override;

private:
	/* Reads levels and their quantum numbers from CHIANTI, as well as the A coefficients
	 * between them. The listed levels are only those involved in spontaneous transitions, up to
	 * n = 5: 1s 2s 2p 2p 3s 3p 3p 3d 3d 4s 4p 4p 4d 4d 4f 4f 5s 5p 5p 5d 5d 5f 5f 5g 5g, and
	 * are also j-resolved (Hence the duplicates shown in this sentence). */
	void readCHIANTI();

	/* Returns energy of a level read in from CHIANTI, given the principal (n), angular momentum
	 * (l) and total momentum (j) quantum numbers. */
	double energy_CHIANTI(int n, int l, int j);

	/* Return the A coefficient between two levels read in from the CHIANTI database. */
	double einsteinAijCHIANTI(int ni, int li, int ji, int nf, int lf, int jf);

	/* Reads collisional data from Anderson+2002 (J. Phys. B: At., Mol. Opt. Phys., 2002, 35,
	 * 1613) */
	void readCollisionData();
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
