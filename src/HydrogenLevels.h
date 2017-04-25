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
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
