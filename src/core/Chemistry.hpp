#ifndef CORE_CHEMISTRY_HPP
#define CORE_CHEMISTRY_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

#include <map>
#include <vector>

class Chemistry
{
protected:
	void addReaction(const std::string& reactionName,
	                 const std::vector<std::string>& reactantNamev,
	                 const Array& reactantStoichv,
	                 const std::vector<std::string>& productNamev,
	                 const Array& productStoichv);

	void addConserved(const std::vector<std::string>& speciesNamev,
	                  const Array& coefficientv);

public:
	virtual ~Chemistry() = default;

	int reactionIndex(const std::string& reactionName) const;
	int numReactions() const { return _reactionv.size(); }
	int numConserved() const { return _conservationv.size(); }
	int numSpecies() const { return _numSpecies; }

	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;
	EVector evaluateFv(const EVector& nv, const EVector& rateCoeffv) const;
	EMatrix evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const;
	EMatrix conservEqvv() const { return _conservEqvv; }

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

	EMatrix _rStoichvv, _netStoichvv, _conservEqvv;
	size_t _numSpecies;
	int _numReactions, _numConserved;

	typedef struct Reaction
	{
		Reaction(const EVector& rv, const EVector& pv) : _rv(rv), _pv(pv) {}
		EVector _rv, _pv;
	} Reaction;
	std::vector<Reaction> _reactionv;
	std::vector<EVector> _conservationv;
	std::map<std::string, int> _reactionIndexm;
};

#endif // CORE_CHEMISTRY_HPP
