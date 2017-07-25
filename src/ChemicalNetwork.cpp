#include "ChemicalNetwork.h"
#include "Error.h"
#include "IonizationBalance.h"

#include <iostream>

namespace
{
const int numSpecies = 4;
const int numConserved = 2;
}

// TODO: need general system for the whole code
const std::map<std::string, int> ChemicalNetwork::speciesIndex = {
                {"e-", 0}, {"H+", 1}, {"H", 2}, {"H2", 3}};

ChemicalNetwork::ChemicalNetwork()
{
	// Add reactions like this:
	// _reactionv.emplace_back(leftside{ne, np, nH, nH2}, rightside{ne, np, nH, nH2})

	// Photoionization
	// H + gamma -> e- + H+
	//_reactionv.emplace_back(Array{0, 0, 1, 0}, Array{1, 1, 0, 0});
	addReaction({"H"}, {1}, {"e-", "H+"}, {1, 1});

	// TODO: add collisional ionization

	// Radiative recombination
	// e- + H+ -> H + gamma
	//_reactionv.emplace_back(Array{1, 1, 0, 0}, Array{0, 0, 1, 0});
	addReaction({"e-", "H+"}, {1, 1}, {"H"}, {1});

	// Dissociation after excitation
	// H2 -> H + H
	//_reactionv.emplace_back(Array{0, 0, 0, 1}, Array{0, 0, 2, 0});
	addReaction({"H2"}, {1}, {"H"}, {2});

	_numReactions = _reactionv.size();
}

void ChemicalNetwork::addReaction(const std::vector<std::string>& reactantNamev,
                                  const Array& reactantStoichv,
                                  const std::vector<std::string>& productNamev,
                                  const Array& productStoichv)
{
	if (reactantNamev.size() != reactantStoichv.size())
		Error::runtime("Error adding reaction: reactantNamev and reactantStoichv size "
		               "mismatch");
	if (productNamev.size() != productStoichv.size())
		Error::runtime("Error adding reaction: productNamev and productStoichv size "
		               "mismatch");

	Array leftSidev(numSpecies);
	for (size_t r = 0; r < reactantNamev.size(); r++)
		leftSidev[speciesIndex.at(reactantNamev[r])] = reactantStoichv[r];

	Array rightSidev(numSpecies);
	for (size_t p = 0; p < productNamev.size(); p++)
		rightSidev[speciesIndex.at(productNamev[p])] = productStoichv[p];

	_reactionv.emplace_back(leftSidev, rightSidev);
}

EMatrix ChemicalNetwork::reactantStoichvv() const
{
	EMatrix r(numSpecies, _numReactions);
	for (int j = 0; j < _numReactions; j++)
		r.col(j) = _reactionv[j]._rv;
	return r;
}

EMatrix ChemicalNetwork::productStoichvv() const
{
	EMatrix p(numSpecies, _numReactions);
	for (int j = 0; j < _numReactions; j++)
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
	EVector k(_numReactions);

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
