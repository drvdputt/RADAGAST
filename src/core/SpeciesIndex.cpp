#include "SpeciesIndex.hpp"

const std::vector<std::string> SpeciesIndex::common4{"e-", "H+", "H", "H2"};

SpeciesIndex::SpeciesIndex(const std::vector<std::string>& namev) : _namev{namev}
{
	for (const std::string& s : namev)
		_indexMap.insert({s, _indexMap.size()});
}

void SpeciesIndex::addSpecies(const std::string& name)
{
	_namev.emplace_back(name);
	_indexMap.insert({name, _indexMap.size()});
}

int SpeciesIndex::index(const std::string& name) const
{
	auto it = _indexMap.find(name);
	if (it != _indexMap.end())
		return it->second;
	else
		return -1;
}

EVector SpeciesIndex::linearCombination(const std::vector<std::string>& namev,
                                        const Array& weightv) const
{
	EVector v = EVector::Zero(size());
	for (size_t i = 0; i < namev.size(); i++)
		v(index(namev[i])) += weightv[i];
	return v;
}

SpeciesVector::SpeciesVector(const SpeciesIndex& speciesIndex)
                : _index{&speciesIndex}, _ine{speciesIndex.index("e-")},
                  _inp{speciesIndex.index("H+")}, _inH{speciesIndex.index("H")},
                  _inH2{speciesIndex.index("H2")}, _nv{EVector::Zero(speciesIndex.size())}
{
}
