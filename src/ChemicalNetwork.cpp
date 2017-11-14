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
	addReaction("H photoionization", {"H"}, {1}, {"e-", "H+"}, {1, 1});

	// Ionization by collision with electron
	// H + e- -> H+ + 2e-
	addReaction("H collisional ionization", {"H", "e-"}, {1, 1}, {"H+", "e-"}, {1, 2});

	// Radiative recombination
	// e- + H+ -> H + gamma
	addReaction("H radiative recombination", {"e-", "H+"}, {1, 1}, {"H"}, {1});

	// Dissociation after excitation
	// H2 -> H + H
	addReaction("H2 dissociation", {"H2"}, {1}, {"H"}, {2});

	/* H2 formation on grain surfaces H + H -> H2 Need to put 1 here for the H stoichiometry,
	   because the rate scales with nH instead of nH^2 (see rateCoeffv()). By then using twice
	   (not the word 'half' in the label) the reaction rate, we end up with an equal tempo of H2
	   formation, but on that scales only linearly with nH. */
	addReaction("half H2 formation", {"H"}, {1}, {"H2"}, {.5});

	_numReactions = _reactionv.size();
}

EVector ChemicalNetwork::rateCoeffv(double T, const Array& frequencyv,
                                    const Array& specificIntensityv, double kDissFromH2Levels,
                                    double kH2FormationGrain) const
{
	EVector k(_numReactions);
	k(reactionIndex("H photoionization")) =
	                Ionization::photoRateCoeff(frequencyv, specificIntensityv);
	k(reactionIndex("H collisional ionization")) = Ionization::collisionalRateCoeff(T);
	k(reactionIndex("H radiative recombination")) = Ionization::recombinationRateCoeff(T);
	k(reactionIndex("H2 dissociation")) = kDissFromH2Levels;
	/* The rate of this reaction is twice the H2 formation rate. See comment in constructor
	   implementation. */
	k(reactionIndex("half H2 formation")) = 2 * kH2FormationGrain;
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
