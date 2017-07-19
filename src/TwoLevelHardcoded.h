#ifndef GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_
#define GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_

#include "LevelDataProvider.h"

class TwoLevelHardcoded : public LevelDataProvider
{
public:
	TwoLevelHardcoded();

	int numLv() const override;
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(double T, double ne, double np) const override;
	EVector sourcev(double T, double ne, double np) const override;
	EVector sinkv(double T, double ne, double np) const override;
private:
	Eigen::Vector2d the_ev;
	Eigen::Vector2d the_gv;
};

#endif /* GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_ */
