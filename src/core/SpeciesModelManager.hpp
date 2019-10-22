#ifndef CORE_SPECIESMODELMANAGER_HPP
#define CORE_SPECIESMODELMANAGER_HPP

#include "H2Data.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"

class HData;
class H2Data;

class SpeciesModelManager
{
public:
	SpeciesModelManager(const std::string& hOption, const std::string& h2Option);
	std::unique_ptr<HModel> makeHModel() const;
	std::unique_ptr<H2Model> makeH2Model(/* heuristic parameters */) const;

private:
	std::unique_ptr<HData> _hData;
	std::unique_ptr<H2Data> _h2Data;
};

#endif // CORE_SPECIESMODELMANAGER_HPP
