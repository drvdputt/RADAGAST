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
	EVector sinkv(double T, double n, double ne, double np) const override;
private:
	EVector the_ev{EVector::Zero(2)};
	EVector the_gv{EVector::Zero(2)};
};

#endif /* GASMODULE_GIT_SRC_TWOLEVELHARDCODED_H_ */
