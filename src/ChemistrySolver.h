#ifndef GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_
#define GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_

#include "EigenAliases.h"

#include <memory>

class ChemicalNetwork;

class ChemistrySolver
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
	ChemistrySolver(const EMatrix& reactantStoichvv, const EMatrix& productStoichvv,
	                const EMatrix& conservationCoeffvv, const EVector& qvv);

	/** A constructor which relies on a ChemicalNetwork object to provide a consistent set of data.
	    The chemical network is in a separate implementation because it will take an unknown amount
	    of arguments to get the rate coefficients. */
	ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn);

	~ChemistrySolver();

	/** Solves the chemical network given a certain rate coefficient vector (indexed on the
	    reactions). An initial value can be given. A vector containing updated densities is
	    returned. The best way to make sure that the rate coefficients are compatible */
	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;

	const ChemicalNetwork* chemicalNetwork() const { return _cn.get(); }

private:
	// DESIGN NOTE:
	/** The calculation of the rate coefficients can be delegated to another object, as the arguments
	    required to calculate them could be anything from the rest of the simulation. This object
	    will provide data that is self-consistent (same number of indices, and in the same order,
	    etc.) To get a set of rate coefficients which is consitent with the data given at setup, the
	    user can then simply ask for the ChemicalNetwork and call its rateCoeffv(<many arguments>) function.
	    One can also opt to not use a chemicalnetwork object, and just fill in everything in
	    the constructor manually. The fact that one then needs to remember what was originally filled
	    in when providing the system with rate coefficients, illustrates why it is useful to have
	    this extra object. */
	std::unique_ptr<const ChemicalNetwork> _cn;
	EMatrix _rvv, _pvv, _cvv;
};

#endif /* GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_ */
