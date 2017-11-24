#include "SpeciesIndex.h"

const std::map<std::string, int> SpeciesIndex::indexMap =
                createSpeciesIndexm({"e-", "H+", "H", "H2"});

std::map<std::string, int> SpeciesIndex::createSpeciesIndexm(const std::vector<std::string>& namev)
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

int SpeciesIndex::index(const std::string& s) { return indexMap.at(s); }

int SpeciesIndex::ine() { return indexMap.at("e-"); }

int SpeciesIndex::inp() { return indexMap.at("H+"); }

size_t SpeciesIndex::size() { return indexMap.size(); }
