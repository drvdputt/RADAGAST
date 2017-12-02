#include "SpeciesIndex.h"

const std::map<std::string, int> SpeciesIndex::_indexMap =
                createSpeciesIndexm({"e-", "H+", "H", "H2"});

std::map<std::string, int>
SpeciesIndex::createSpeciesIndexm(const std::vector<std::string>& namev)
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

int SpeciesIndex::index(const std::string& s) { return _indexMap.at(s); }

int SpeciesIndex::ine() { return _indexMap.at("e-"); }

int SpeciesIndex::inp() { return _indexMap.at("H+"); }

int SpeciesIndex::inH() { return _indexMap.at("H"); }

int SpeciesIndex::inH2() { return _indexMap.at("H2"); }

size_t SpeciesIndex::size() { return _indexMap.size(); }
