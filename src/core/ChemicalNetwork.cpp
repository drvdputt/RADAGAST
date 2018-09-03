#include "ChemicalNetwork.h"
#include "SpeciesIndex.h"

#include <iostream>

ChemicalNetwork::ChemicalNetwork() : _numSpecies{SpeciesIndex::size()} {}

void ChemicalNetwork::addReaction(const std::string& reactionName,
                                  const std::vector<std::string>& reactantNamev,
                                  const Array& reactantStoichv,
                                  const std::vector<std::string>& productNamev,
                                  const Array& productStoichv)
{
	// Give the reaction a number, and put its name in the map
	_reactionIndexm.emplace(reactionName, _reactionIndexm.size());
	_reactionv.emplace_back(
	                SpeciesIndex::makeFullCoefficientv(reactantNamev, reactantStoichv),
	                SpeciesIndex::makeFullCoefficientv(productNamev, productStoichv));
}

void ChemicalNetwork::addConserved(const std::vector<std::string>& speciesNamev,
                                   const Array& coefficientv)
{
	_conservationv.emplace_back(
	                SpeciesIndex::makeFullCoefficientv(speciesNamev, coefficientv));
}

int ChemicalNetwork::reactionIndex(const std::string& reactionName) const
{
	return _reactionIndexm.at(reactionName);
}

EMatrix ChemicalNetwork::reactantStoichvv() const
{
	EMatrix r(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
	{
		// Each column represents a reaction
		r.col(j) = _reactionv[j]._rv;
	}
	return r;
}

EMatrix ChemicalNetwork::productStoichvv() const
{
	EMatrix p(_numSpecies, _reactionv.size());
	for (size_t j = 0; j < _reactionv.size(); j++)
	{
		// Each column represents a reaction
		p.col(j) = _reactionv[j]._pv;
	}
	return p;
}

EMatrix ChemicalNetwork::conservationCoeffvv() const
{
	EMatrix c(numConserved(), numSpecies());
	for (size_t q = 0; q < _conservationv.size(); q++)
	{
		// Each row represents a conserved quantity
		c.row(q) = _conservationv[q];
	}
	return c;
}
