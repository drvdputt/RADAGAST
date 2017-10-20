#ifndef H2FROMFILES_H_
#define H2FROMFILES_H_

#include "LevelDataProvider.h"

class H2FromFiles : public LevelDataProvider
{
public:
	H2FromFiles();
	~H2FromFiles();
	int numLv() const override;
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(double T, double ne, double np) const override;
	EVector sourcev(double T, double ne, double np) const override;
	EVector sinkv(double T, double n, double ne, double np) const override;
};

#endif /* H2FROMFILES_H_ */
