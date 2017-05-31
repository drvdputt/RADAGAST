#ifndef _TWOLEVEL_H_
#define _TWOLEVEL_H_

#include "NLevel.h"

#include <Eigen/Dense>

class TwoLevel : public NLevel
{
public:
	/* Creates an object that represents a two-level (component of a) medium. The level
	 population equilibrium will be calculated using the frequency grid supplied as argument of
	 the constructor. */
	TwoLevel(const Array& frequencyv);

protected:
	int makeNLv() const override;
	Eigen::VectorXd makeEv() const override;
	Eigen::VectorXd makeGv() const override;
	Eigen::MatrixXd makeAvv() const override;

	Eigen::MatrixXd prepareCollisionMatrix(double T, double electronDensity,
	                                       double protonDensity) const override;
};

#endif /* _TWOLEVEL_H_ */
