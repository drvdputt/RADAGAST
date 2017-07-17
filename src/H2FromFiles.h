#ifndef H2FROMFILES_H_
#define H2FROMFILES_H_

#include "LevelDataProvider.h"

class H2FromFiles: public LevelDataProvider
{
public:
	H2FromFiles();
	~H2FromFiles();
	int numLv() const override;
	Eigen::VectorXd ev() const override;
	Eigen::VectorXd gv() const override;
	Eigen::MatrixXd avv() const override;
	Eigen::MatrixXd extraAvv() const override;
	Eigen::MatrixXd cvv(double T, double ne, double np) const override;
	Eigen::VectorXd sourcev(double T, double ne, double np) const override;
	Eigen::VectorXd sinkv(double T, double ne, double np) const override;
};

#endif /* H2FROMFILES_H_ */
