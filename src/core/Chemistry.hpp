#ifndef CORE_CHEMISTRY_HPP
#define CORE_CHEMISTRY_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

#include <map>
#include <vector>



class Chemistry
{
public:
	virtual ~Chemistry() = default;
	void addReaction(const std::string& reactionName,
	                 const std::vector<std::string>& reactantNamev,
	                 const Array& reactantStoichv,
	                 const std::vector<std::string>& productNamev,
	                 const Array& productStoichv);

	/** This should be called after reactions have been added */
	void prepareCoefficients();

	int reactionIndex(const std::string& reactionName) const;
	int numReactions() const { return _reactionv.size(); }
	int numSpecies() const { return _numSpecies; }

	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;
	EVector evaluateFv(const EVector& nv, const EVector& rateCoeffv) const;
	EMatrix evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const;

private:
	// help setting the members
	EMatrix makeReactantStoichvv() const;
	EMatrix makeProductStoichvv() const;
	EMatrix makeConservationCoeffvv() const;

	EVector solveTimeDep(const EVector& rateCoeffv, const EVector& n0v) const;
	double densityProduct(const EVector& nv, size_t r) const;
	double densityProductDerivative(const EVector& nv, int r, int j) const;

	// I'm currently combining the members from ChemistrySolver and ChemicalNetwork. Should
	// refactor out the redundancies later, after tests have been written/ported.

	typedef struct Reaction
	{
		Reaction(const EVector& rv, const EVector& pv) : _rv(rv), _pv(pv) {}
		EVector _rv, _pv;
	} Reaction;
	std::vector<Reaction> _reactionv;
	std::map<std::string, int> _reactionIndexm;

	EMatrix _rStoichvv, _netStoichvv;
	size_t _numSpecies;
	int _numReactions;
};

#endif // CORE_CHEMISTRY_HPP
