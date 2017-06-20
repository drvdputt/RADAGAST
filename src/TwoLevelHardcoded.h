#ifndef GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_
#define GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_

#include "LevelDataProvider.h"

class TwoLevelHardcoded : public LevelDataProvider
{
public:
	TwoLevelHardcoded();

	int numLv() const override;
	Eigen::VectorXd ev() const override;
	Eigen::VectorXd gv() const override;
	Eigen::MatrixXd avv() const override;
	Eigen::MatrixXd extraAvv() const override;

	Eigen::MatrixXd cvv(double T, double ne, double np) const override;
	Eigen::VectorXd alphav(double T) const override;
private:
	Eigen::Vector2d the_ev;
	Eigen::Vector2d the_gv;
};

#endif /* GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_ */
