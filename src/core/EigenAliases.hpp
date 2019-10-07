#ifndef CORE_EIGENALIASES_H_
#define CORE_EIGENALIASES_H_

#include <Eigen/Dense>

using EMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using EMatrix_bool = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
// Row major version can be useful for interfacing with e.g. GSL
using EMatrixRM = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using EVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EArray = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

#endif /* CORE_EIGENALIASES_H_ */