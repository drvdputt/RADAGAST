#ifndef GASMODULE_GIT_SRC_EIGENALIASES_H_
#define GASMODULE_GIT_SRC_EIGENALIASES_H_

#include <Eigen/Dense>

using EMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using EVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using EArray = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Derived>
void compareMatrices(const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<Derived>& b,
                     double tolerance)
{
	// Calculate relative difference element-wise
	EMatrix relDiff = (a - b).array() / b.array();
	for (int i = 0; i < relDiff.size(); i++)
	{
		// Take out the nans (due to divide by zero)
		auto* pointer = relDiff.data() + i;
		auto value = *pointer;
		*pointer = isnan(value) ? 0 : value;
	}
	assert((relDiff.cwiseAbs().array() < tolerance).all());
}

#endif /* GASMODULE_GIT_SRC_EIGENALIASES_H_ */
