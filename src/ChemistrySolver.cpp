#include "ChemistrySolver.h"
#include "ChemicalNetwork.h"

ChemistrySolver::ChemistrySolver(const EMatrix& reactantStoichvv, const EMatrix& productStoichvv,
                                 const EMatrix& conservationCoeffvv, const EVector& qv)
                : _rvv(reactantStoichvv), _pvv(productStoichvv), _cvv(conservationCoeffvv)
{
}

ChemistrySolver::ChemistrySolver(std::unique_ptr<const ChemicalNetwork> cn) : _cn(std::move(cn))
{
	_rvv = _cn->reactantStoichvv();
	_pvv = _cn->productStoichvv();
	_cvv = _cn->conservationCoeffvv();
}

ChemistrySolver::~ChemistrySolver() = default;

EVector ChemistrySolver::solveBalance(const EVector& rateCoeffv, const EVector& n0v) const
{
	// TODO: Newton-Raphson
	EVector r(4);
	return r;
}
