#include "SpeciesIndex.hpp"
#include "Error.hpp"

// const std::map<std::string, int> SpeciesIndex::_indexMap =
//                 createSpeciesIndexm({"e-", "H+", "H", "H2"});

// std::map<std::string, int>
// SpeciesIndex::createSpeciesIndexm(const std::vector<std::string>& namev)
// {
// 	std::map<std::string, int> speciesIndex;
// 	int index{0};
// 	for (const auto& name : namev)
// 	{
// 		speciesIndex.insert({name, index});
// 		index++;
// 	}
// 	return speciesIndex;
// }

// int SpeciesIndex::index(const std::string& s) { return _indexMap.at(s); }

// int SpeciesIndex::ine() { return _indexMap.at("e-"); }

// int SpeciesIndex::inp() { return _indexMap.at("H+"); }

// int SpeciesIndex::inH() { return _indexMap.at("H"); }

// int SpeciesIndex::inH2() { return _indexMap.at("H2"); }

// size_t SpeciesIndex::size() { return _indexMap.size(); }

// EVector SpeciesIndex::makeFullCoefficientv(const std::vector<std::string>& namev,
//                                            const Array& coefficientv)
// {
// 	Error::equalCheck("Lengths of list of species names and vector of coefficients",
// 	                  namev.size(), coefficientv.size());

// 	EVector fullCoefficientv = EVector::Zero(size());
// 	for (size_t r = 0; r < namev.size(); r++)
// 		fullCoefficientv(index(namev[r])) += coefficientv[r];

// 	return fullCoefficientv;
// }

NewSpeciesIndex::NewSpeciesIndex(const std::vector<std::string>& namev) : _namev{namev}
{
	for (const std::string& s : namev)
		_indexMap.insert({s, _indexMap.size()});
}

void NewSpeciesIndex::addSpecies(const std::string& name)
{
	_namev.emplace_back(name);
	_indexMap.insert({name, _indexMap.size()});
}

int NewSpeciesIndex::index(const std::string& name) const
{
	auto it = _indexMap.find(name);
	if (it != _indexMap.end())
		return it->second;
	else
		return -1;
}

EVector NewSpeciesIndex::unitVector(const std::string& name)
{
	EVector v = EVector::Zero(size());
	v(index(name)) = 1;
	return v;
}

SpeciesVector::SpeciesVector(const NewSpeciesIndex& speciesIndex)
                : _index{&speciesIndex}, _ine{speciesIndex.index("e-")},
                  _inp{speciesIndex.index("H+")}, _inH{speciesIndex.index("H")},
                  _inH2{speciesIndex.index("H2")}, _nv{EVector::Zero(speciesIndex.size())}
{
}
