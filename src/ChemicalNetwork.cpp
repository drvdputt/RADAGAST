#include "ChemicalNetwork.h"
#include "IonizationBalance.h"

namespace
{
const int numReactions = 4;
const int numSpecies = 4;
const int numConserved = 2;
// TODO temp hack, need general system for the whole code
//int ine = 0, inp = 1, inH = 2, inH2 = 3;
// Species are (electron, proton, atomic H, molecular H2)
}

ChemicalNetwork::ChemicalNetwork()
{
	// Add reactions like this:
	// _reactionv.emplace_back(leftside{ne, np, nH, nH2}, rightside{ne, np, nH, nH2})

	// Photoionization
	// H + gamma -> ne + np
	_reactionv.emplace_back(Array{0, 0, 1, 0}, Array{1, 1, 0, 0});

	// TODO: add collisional ionization

	// Radiative recombination
	// ne + np -> H + gamma
	_reactionv.emplace_back(Array{1, 1, 0, 0}, Array{0, 0, 1, 0});

	// Dissociation after excitation
	// H2 -> H + H
	_reactionv.emplace_back(Array{0, 0, 0, 1}, Array{0, 0, 2, 0});
}

EMatrix ChemicalNetwork::reactantStoichvv() const
{
	EMatrix r(numSpecies, numReactions);
	for (int j = 0; j < numReactions; j++)
		r.col(j) = _reactionv[j]._rv;
	return r;
}

EMatrix ChemicalNetwork::productStoichvv() const
{
	EMatrix p(numSpecies, numReactions);
	for (int j = 0; j < numReactions; j++)
		p.col(j) = _reactionv[j]._pv;
	return p;
}

EMatrix ChemicalNetwork::conservationCoeffvv() const
{
	EMatrix c(numConserved, numSpecies);

	// Conservation of H nuclei (single protons)
	EVector protonEq(numSpecies);
	protonEq << 0, 1, 1, 2;
	c.row(0) = protonEq;

	// Conservation of electrons
	EVector electronEq(numSpecies);
	electronEq << 1, 0, 1, 2;
	c.row(1) = electronEq;

	return c;
}

EVector ChemicalNetwork::rateCoeffv(double T, const Array& frequencyv,
                                    const Array& specificIntensityv, double kFromH2Levels) const
{
	EVector k(numReactions);

	// Photoionization
	// H + gamma -> ne + np
	k(0) = Ionization::photoRateCoeff(frequencyv, specificIntensityv);

	// TODO: add collisional ionization

	// Radiative recombination
	// ne + np -> H + gamma
	k(1) = Ionization::recombinationRateCoeff(T);

	// Dissociation after excitation
	// H2 -> H + H
	k(2) = kFromH2Levels;

	return k;
}
