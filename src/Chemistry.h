#ifndef GASMODULE_GIT_SRC_CHEMISTRY_H_
#define GASMODULE_GIT_SRC_CHEMISTRY_H_

#include <Eigen/Dense>

class Chemistry
{
public:
	/** Create a chemistry network with coefficients specified in the following matrices:
	    reactantStoichvv contains the stoichiometry for each reactant (the left side of the
	    reaction). productStoichvv contains the stoichiometric numbers for the reaction products
	    (the right side). The rows of these matrices are indexed per species, while each column
	    stands for a certain reaction. It will be assumed that the reaction rates scale with the
	    product of the reactant densities (to the power of their stoich. number on the left hand
	    side). This makes it easy to compute the Jacobian of the time derivative of the
	    densities. The last two arguments are a matrix and a vector representing conservation
	    equations. Each row represents the conservation of a certain quantity q;
	    conservationCoeffvv should contain the amount of q contained in each species of the
	    chemical network. The complementary vector qvv indicates the desired value of each q. By
	    multiplying the convervation coefficient matrix with the species density vector, we
	    obtain as such our set of conservation equations. */
	Chemistry(const Eigen::MatrixXd& reactantStoichvv, const Eigen::MatrixXd& productStoichvv,
	          const Eigen::MatrixXd& conservationCoeffvv, const Eigen::VectorXd& qvv);

	/** Solves the chemical network given a certain rate coefficient vector (indexed on the
	    reactions). An initial value can be given. A vector containing updated densities is
	    returned. */
	Eigen::VectorXd solveBalance(const Eigen::VectorXd rateCoeffv, const Eigen::VectorXd n0v);

private:
	Eigen::MatrixXd _rvv, _pvv, _cvv;
	Eigen::VectorXd _qvv;
};

#endif /* GASMODULE_GIT_SRC_CHEMISTRY_H_ */
