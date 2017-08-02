#include "ChemicalNetwork.h"
#include "Error.h"
#include "IonizationBalance.h"

#include <iostream>

namespace
{
const int NUMSPECIES = 4;
const int NUMCONSERVED = 2;
}

const std::map<std::string, int> ChemicalNetwork::speciesIndexm =
                createSpeciesIndexm({"e-", "H+", "H", "H2"});

std::map<std::string, int>
ChemicalNetwork::createSpeciesIndexm(const std::vector<std::string>& namev)
{
	std::map<std::string, int> speciesIndex;
	int index{0};
	for (const auto& name : namev)
	{
		speciesIndex.insert({name, index});
		index++;
	}
	return speciesIndex;
}

ChemicalNetwork::ChemicalNetwork()
{
	// Photoionization
	// H + gamma -> e- + H+
	//_reactionv.emplace_back(Array{0, 0, 1, 0}, Array{1, 1, 0, 0});
	addReaction("H photoionization", {"H"}, {1}, {"e-", "H+"}, {1, 1});

	// Ionization by collision with electron
	// H + e- -> H+ + 2e-
	addReaction("H collisional ionization", {"H", "e-"}, {1, 1}, {"H+", "e-"}, {1, 2});

	// Radiative recombination
	// e- + H+ -> H + gamma
	//_reactionv.emplace_back(Array{1, 1, 0, 0}, Array{0, 0, 1, 0});
	addReaction("H radiative recombination", {"e-", "H+"}, {1, 1}, {"H"}, {1});

	// Dissociation after excitation
	// H2 -> H + H
	//_reactionv.emplace_back(Array{0, 0, 0, 1}, Array{0, 0, 2, 0});
	addReaction("H2 dissociation", {"H2"}, {1}, {"H"}, {2});

	_numReactions = _reactionv.size();
}

EVector ChemicalNetwork::rateCoeffv(double T, const Array& frequencyv,
                                    const Array& specificIntensityv, double kFromH2Levels) const
{
	EVector k(_numReactions);

	// Photoionization
	// H + gamma -> ne + np
	k(reactionIndex("H photoionization")) =
	                Ionization::photoRateCoeff(frequencyv, specificIntensityv);

	// Ionization by collision with electron
	k(reactionIndex("H collisional ionization")) = Ionization::collisionalRateCoeff(T);

	// Radiative recombination
	// ne + np -> H + gamma
	k(reactionIndex("H radiative recombination")) = Ionization::recombinationRateCoeff(T);

	// Dissociation after excitation
	// H2 -> H + H
	k(reactionIndex("H2 dissociation")) = kFromH2Levels;

	return k;
}

void ChemicalNetwork::addReaction(const std::string& reactionName,
                                  const std::vector<std::string>& reactantNamev,
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

	_reactionIndexm.emplace(reactionName, _reactionIndexm.size());

	EVector leftSidev = EVector::Zero(NUMSPECIES);
	for (size_t r = 0; r < reactantNamev.size(); r++)
		leftSidev(speciesIndexm.at(reactantNamev[r])) += reactantStoichv[r];

	EVector rightSidev = EVector::Zero(NUMSPECIES);
	for (size_t p = 0; p < productNamev.size(); p++)
		rightSidev(speciesIndexm.at(productNamev[p])) += productStoichv[p];

	_reactionv.emplace_back(leftSidev, rightSidev);
}

int ChemicalNetwork::reactionIndex(const std::string& reactionName) const
{
	return _reactionIndexm.at(reactionName);
}

EMatrix ChemicalNetwork::reactantStoichvv() const
{
	EMatrix r(NUMSPECIES, _numReactions);
	for (int j = 0; j < _numReactions; j++)
		r.col(j) = _reactionv[j]._rv;
	return r;
}

EMatrix ChemicalNetwork::productStoichvv() const
{
	EMatrix p(NUMSPECIES, _numReactions);
	for (int j = 0; j < _numReactions; j++)
		p.col(j) = _reactionv[j]._pv;
	return p;
}

EMatrix ChemicalNetwork::conservationCoeffvv() const
{
	EMatrix c(NUMCONSERVED, NUMSPECIES);

	// Some index which should just go from 0 to NUMCONSERVED - 1
	int index{0};

	// Commonly used
	int iH = speciesIndexm.at("H");
	int iH2 = speciesIndexm.at("H2");

	// Conservation of H nuclei (single protons)
	EVector protonEq{EVector::Zero(NUMSPECIES)};
	protonEq(speciesIndexm.at("H+")) = 1;
	protonEq(iH) = 1;
	protonEq(iH2) = 2;
	c.row(index) = protonEq;
	index++;

	// Conservation of electrons
	EVector electronEq{EVector::Zero(NUMSPECIES)};
	electronEq(speciesIndexm.at("e-")) = 1;
	electronEq(iH) = 1;
	electronEq(iH2) = 2;
	c.row(index) = electronEq;
	index++;

	return c;
}
