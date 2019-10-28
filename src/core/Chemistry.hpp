#ifndef CORE_CHEMISTRY_HPP
#define CORE_CHEMISTRY_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "SpeciesIndex.hpp"

#include <map>
#include <vector>

/** This class can store a list of reactions, and solve the chemical equilibrium based on the
    coefficients of the reactions and their reaction rates. */
class Chemistry
{
public:
	virtual ~Chemistry() = default;

	/** TODO & CAUTION: Do not call after adding reactions! Things WILL break. TODO: change
	    the `Reaction' struct to store the reaction information in a way that doesn't rely
	    on the size of the system (just names and stoich should be fine, prepareCoefficients
	    can then convert everything properly). Do this after the setup and storage and
	    passing around of the new SpeciesIndex has been implemented. */
	void registerSpecies(const std::vector<std::string>& namev);

	/** Read-only acces to a SpeciesIndex object, in case detailed information about the
	    SpeciesVector is needed. Useful for making initial guesses. */
	const SpeciesIndex& speciesIndex() const { return _speciesIndex; }

	/** Function to provide a clear syntax for adding reactions in the setup of the chemical
	    network. Each reaction is given a number, and the reaction is added to the reaction
	    index using the given name as a key. Subclasses typically implement a constructor
	    which calls addReaction and prepareCoefficients. Manually adding an extra reaction
	    is still possible afterwards, but make sure to call prepareCoefficients again. On
	    the other hand, I disallowed adding extra species, because that will break the way
	    the current reactions are stored (as vectors in the species space)..*/
	void addReaction(const std::string& reactionName,
	                 const std::vector<std::string>& reactantNamev,
	                 const Array& reactantStoichv,
	                 const std::vector<std::string>& productNamev,
	                 const Array& productStoichv);

	/** This should be called after reactions have been added */
	void prepareCoefficients();

	/** Look up the index of a reaction, based on the name that was given in to addReaction. */
	int reactionIndex(const std::string& reactionName) const;

	/** The number of reactions. This determines the size of the reaction rate vector */
	int numReactions() const { return _reactionv.size(); }

	/** The number of species. This determines the size of the species density vector */
	int numSpecies() const { return _numSpecies; }

	/** Solves the chemical network given a certain rate coefficient vector (indexed on the
	    reactions). An initial value n0v can be given. A vector containing updated densities
	    is returned. */
	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;

	/** Evaluate the rate of change for each species [cm-3 s-1]. */
	EVector evaluateFv(const EVector& nv, const EVector& rateCoeffv) const;

	/** Evaluate the Jacobian of Fv [s-1]. Every column j is the derivative of Fv towards
	    n_j. */
	EMatrix evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const;

private:
	// Turn list of reactions into coefficient matrices
	EMatrix makeReactantStoichvv() const;
	EMatrix makeProductStoichvv() const;

	/** Solve the chemistry by evolving the system until equilibrium. */
	EVector solveTimeDep(const EVector& rateCoeffv, const EVector& n0v) const;

	/** Calculate the density factor needed to calculate the speed of reaction r. Formula:
	    Product_i n_i ^ Rir, where n_i are the elements of nv, and Rir = _rStoichvv(i,
	    r). */
	double densityProduct(const EVector& nv, size_t r) const;

	/** Calculate the derivative of the density product for reaction r with respect to the
	    density j. Formula: (Rjr - 1) * n_j^{Rjr - 1} * Product_{i != j} n_i ^ Rir */
	double densityProductDerivative(const EVector& nv, int r, int j) const;

	// Keep track of index for each species name.
	SpeciesIndex _speciesIndex;

	typedef struct Reaction
	{
		Reaction(const EVector& rv, const EVector& pv) : _rv(rv), _pv(pv) {}
		EVector _rv, _pv;
	} Reaction;
	std::vector<Reaction> _reactionv;
	std::map<std::string, int> _reactionIndexm;

	// Filled in by prepareCoefficients()
	EMatrix _rStoichvv, _netStoichvv;
	size_t _numSpecies;
	int _numReactions;
};

#endif // CORE_CHEMISTRY_HPP
