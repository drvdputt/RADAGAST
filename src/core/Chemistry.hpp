#ifndef CORE_CHEMISTRY_HPP
#define CORE_CHEMISTRY_HPP

#include "EigenAliases.hpp"

class Chemistry
{
public:
	virtual ~Chemistry() = default;

	EVector solveBalance(const EVector& rateCoeffv, const EVector& n0v) const;
	EVector evaluateFv(const EVector& nv, const EVector& rateCoeffv) const;
	EMatrix evaluateJvv(const EVector& nv, const EVector& rateCoeffv) const;
	EMatrix conservEqvv() const { return _conservEqvv; }

private:
	EVector solveTimeDep(const EVector& rateCoeffv, const EVector& n0v) const;
	double densityProduct(const EVector& nv, size_t r) const;
	double densityProductDerivative(const EVector& nv, int r, int j) const;

	EMatrix _rStoichvv, _netStoichvv, _conservEqvv;
	size_t _numSpecies;
	int _numReactions, _numConserved;
};

#endif // CORE_CHEMISTRY_HPP
