#ifndef GASMODULE_GIT_SRC_CHEMICALNETWORK_H_
#define GASMODULE_GIT_SRC_CHEMICALNETWORK_H_

#include "EigenAliases.h"
#include "Spectrum.h"

#include <map>
#include <string>
#include <vector>

/** An object of this class will provide a self-consistent set of data that can be handed to a
    ChemistrySolver. Right now, this class contains a concrete implementation prototype, but
    might at some point be subclassed with names such as GloverChemicalNetwork,
    SimpleChemicalNetwork. */
class ChemicalNetwork
{
public:
	/** Sets up a chemical network with a fixed set of reactions, added in the
	    constructor. */
	ChemicalNetwork();

	/** Calculates the rate coefficients for each reaction. They still have to be multiplied
	    with the densities of the reaction products, so watch out. The unit varies, but
	    multiplying with the correct densities (for example *np*ne for radiative
	    recombination) will amount to cm-3 s-1 as expected. Some of the arguments are rates
	    which have been calculated somewhere else. This function simply fills them in into
	    the right spot of the k-vector. */
	EVector rateCoeffv(double T, const Spectrum& specificIntensity,
	                   double kDissFromH2Levels, double kH2FormationGrain) const;

private:
	/** Function to provide a clear syntax for adding reactions in the setup of the chemical
	    network. Each reaction is given a number, and the reaction is added to the reaction
	    index using the given name as a key. */
	void addReaction(const std::string& reactionName,
	                 const std::vector<std::string>& reactantNamev,
	                 const Array& reactantStoichv,
	                 const std::vector<std::string>& productNamev,
	                 const Array& productStoichv);

public:
	/** Look up the index of a reaction. You'll need to look into the source code for the
	    names though... */
	int reactionIndex(const std::string& reactionName) const;

	/** Provides a matrix containing the stoichiometry of each included species s (rows) on
	    the left hand side of each reaction r (columns) of the network. Multiplying this
	    matrix with a column vector of reaction speeds (in cm-3 s-1) will yield the
	    destruction rate of each species. */
	EMatrix reactantStoichvv() const;

	/** Provides a matrix containing the stoichiometry of each included species s (rows) on
	    the right hand side of each reaction r (columns) of the network. Multiplying this
	    matrix with a column vector of reaction speeds (in cm-3 s-1) will yield the creation
	    rate of each species. */
	EMatrix productStoichvv() const;

	/** A matrix which specifies a set of conserved quantities. Each row represents one
	    conservation equation as a set of (positive integer) coefficients, and basically
	    counts the conserved quantity when multiplied with a column vector containing the
	    densities of the species. The exact values of the conserved quantities will be
	    calculated in this way from the initial guess of the chemical abundances. */
	EMatrix conservationCoeffvv() const;

private:
	size_t _numSpecies;

	/** Map reaction names to indices. This is handy if one wants to know what the rate
	    coefficients mean. Reactions can be dynamically added by using the addReaction
	    function. */
	std::map<std::string, int> _reactionIndexm;

	/** Put the reactions in a list for easy processing. The rate coefficients will be
	    calculated individually though (i.e. not in a loop), while respecting the same
	    order. */
	typedef struct Reaction
	{
		Reaction(const EVector& rv, const EVector& pv) : _rv(rv), _pv(pv) {}
		EVector _rv, _pv;
	} Reaction;

	std::vector<Reaction> _reactionv;
	int _numReactions;
};

#endif /* GASMODULE_GIT_SRC_CHEMICALNETWORK_H_ */
