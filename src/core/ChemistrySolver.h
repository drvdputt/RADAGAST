#ifndef GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_
#define GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_

#include "EigenAliases.h"

#include <memory>

class ChemicalNetwork;

/** This class contains an implementation to solve a chemical network. The details of the
    chemical network are provided by the instance of a \c ChemicalNetwork given to the
    constructor of this class. A pointer to the given \c ChemicalNetwork object will be stored.
    There is a getter for this pointer so the client can access it to generate reaction rate
    coefficients. Only the chemical network object knows what the reactions really mean; In
    fact, this class does not even know which element has which index in the solution vector.
    This class just solves the chemical equations oblivious of any underlying details. */
class ChemistrySolver
{
public:
	/** This constructor is pretty difficult to use, but this description nicely describes
	    the responsibilities of the \c ChemicalNetwork class. Create a chemistry network
	    with coefficients specified in the following matrices: reactantStoichvv contains the
	    stoichiometry for each reactant (the left side of the reaction). productStoichvv
	    contains the stoichiometric numbers for the reaction products (the right side). The
	    rows of these matrices are indexed per species, while each column stands for a
	    certain reaction. It will be assumed that the reaction rates scale with the product
	    of the reactant densities (to the power of their stoich. number on the left hand
	    side). This makes it easy to compute the Jacobian of the time derivative of the
	    densities. The last two arguments are a matrix and a vector representing
	    conservation equations. Each row represents the conservation of a certain quantity
	    q; conservationCoeffvv should contain the amount of q contained in each species of
	    the chemical network. The complementary vector qvv indicates the desired value of
	    each q. By multiplying the convervation coefficient matrix with the species density
	    vector, we obtain as such our set of conservation equations. */
	ChemistrySolver(const EMatrix& reactantStoichvv, const EMatrix& productStoichvv,
	                const EMatrix& conservationCoeffvv);

	/** A constructor which relies on a ChemicalNetwork object to provide a consistent set
	    of data. The chemical network is in a separate implementation because it will take
	    an unknown amount of arguments to get the rate coefficients. */
	ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn);

	~ChemistrySolver();

	/** Solves the chemical network given a certain rate coefficient vector (indexed on the
	    reactions). An initial value can be given. A vector containing updated densities is
	    returned. The best way to make sure that the rate coefficients are compatible */
	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;

	/** Returns a pointer to the chemical network that was given at construction. This
	    chemical network can then be used to calculate the rate coefficients in a nicely
	    consistent way, with a set of arguemnts that depends on the type of chemical network
	    used. */
	const ChemicalNetwork* chemicalNetwork() const { return _cn.get(); }

	/** Evaluates the function of which the root needs to be found. This is a combination of
	    net rates being zero, and conservation equations. */
	EVector evaluateFv(const EVector& nv, const EVector& rateCoeffv,
	                   const EVector& conservedQuantityv) const;

	/** Evaluates the jacobian of the function above, i.e. each column is the derivative of
	    Fv with respect to one of the densities. Thanks to the power law prescription of the
	    total rates, this can be calculated analytically. */
	EMatrix evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const;

private:
	/** Calculates the density factor for the reaction r. Multiplying the reaction rate
	    coefficient with this factor gives the total reaction rate. */
	double densityProduct(const EVector& nv, size_t r) const;

	/** Calculates the derivative of the density factor for reaction r, derived to the
	    density of species j. */
	double densityProductDerivative(const EVector& nv, int r, int j) const;

	// DESIGN NOTE:
	/** The calculation of the rate coefficients can be delegated to another object, as the
	    arguments required to calculate them could be anything from the rest of the
	    simulation. This object will provide data that is self-consistent (same number of
	    indices, and in the same order, etc.) To get a set of rate coefficients which is
	    consistent with the data given at setup, the user can then simply ask for the
	    ChemicalNetwork and call its rateCoeffv(<many arguments>) function. One can also opt
	    to not use a chemicalnetwork object, and just fill in everything in the constructor
	    manually. The fact that one then needs to remember what was originally filled in
	    when providing the system with rate coefficients, illustrates why it is useful to
	    have this extra object. */
	std::unique_ptr<const ChemicalNetwork> _cn;
	EMatrix _rStoichvv, _netStoichvv, _conservEqvv;
	size_t _numSpecies;
	int _numReactions, _numConserved;
};

#endif /* GASMODULE_GIT_SRC_CHEMISTRYSOLVER_H_ */