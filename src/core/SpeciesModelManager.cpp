#include "SpeciesModelManager.hpp"
#include "BigH2Model.hpp"
#include "HFromFiles.hpp"
#include "HHardCoded.hpp"
#include "SimpleH2.hpp"

#include <sstream>

SpeciesModelManager::SpeciesModelManager(const std::string& hOption,
                                         const std::string& h2Option)
{
	// Choose hydrogen data
	if (hOption == "hhc")
		_hData = std::make_unique<HydrogenHardcoded>();
	else if (hOption == "hff2")
		_hData = std::make_unique<HFromFiles>(2);
	else if (hOption == "hff4")
		_hData = std::make_unique<HFromFiles>(4);
	else
		_hData = std::make_unique<HFromFiles>();

	// H2 data
	if (h2Option.empty())
		_h2Data = std::make_unique<H2Data>();
	else
	{
		int maxJ, maxV;
		std::istringstream(h2Option) >> maxJ >> maxV;
		if (maxJ < 0 || maxV < 0)
			Error::runtime("moleculeChoice is not of a correct format. It "
			               "should be "
			               "\"maxJ maxV\"");
		_h2Data = std::make_unique<H2Data>(maxJ, maxV);
	}
}

std::unique_ptr<HModel> SpeciesModelManager::makeHModel() const
{
	return std::make_unique<HModel>(_hData.get());
}

std::unique_ptr<H2Model> SpeciesModelManager::makeH2Model(bool simple) const
{
	if (simple)
		return std::make_unique<SimpleH2>();
	else
		return std::make_unique<BigH2Model>(_h2Data.get());
}
