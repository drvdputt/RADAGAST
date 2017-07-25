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

	// Ionization by collision with electron
	// H + e- -> H+ + 2e-
	addReaction({"H", "e-"}, {1, 1}, {"H+", "e-"}, {1, 2});

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

EVector ChemicalNetwork::rateCoeffv(double T, const Array& frequencyv,
                                    const Array& specificIntensityv, double kFromH2Levels) const
{
	EVector k(_numReactions);

	int r = 0;

	// Photoionization
	// H + gamma -> ne + np
	k(r) = Ionization::photoRateCoeff(frequencyv, specificIntensityv);

	// Ionization by collision with electron
	r++;
	k(r) = Ionization::collisionalRateCoeff(T);

	// Radiative recombination
	// ne + np -> H + gamma
	r++;
	k(r) = Ionization::recombinationRateCoeff(T);

	// Dissociation after excitation
	// H2 -> H + H
	r++;
	k(r) = kFromH2Levels;

	return k;
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

	EVector leftSidev = EVector::Zero(numSpecies);
	for (size_t r = 0; r < reactantNamev.size(); r++)
		leftSidev(speciesIndex.at(reactantNamev[r])) += reactantStoichv[r];

	EVector rightSidev= EVector::Zero(numSpecies);
	for (size_t p = 0; p < productNamev.size(); p++)
		rightSidev(speciesIndex.at(productNamev[p])) += productStoichv[p];

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
