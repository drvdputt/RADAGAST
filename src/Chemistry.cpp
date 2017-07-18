#include "Chemistry.h"

Chemistry::Chemistry(const Eigen::MatrixXd& reactantStoichvv,
                     const Eigen::MatrixXd& productStoichvv,
                     const Eigen::MatrixXd& conservationCoeffvv, const Eigen::VectorXd& qvv)
                : _rvv(reactantStoichvv), _pvv(productStoichvv), _cvv(conservationCoeffvv),
                  _qvv(qvv)
{
}
