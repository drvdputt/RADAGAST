#include "SpeciesModelManager.hpp"
#include "HydrogenFromFiles.hpp"
#include "HydrogenHardcoded.hpp"

#include <sstream>

SpeciesModelManager::SpeciesModelManager(const std::string& hOption, const std::string& h2Option)
{
	// Choose hydrogen data
	if (hOption == "hhc")
		_hData = std::make_unique<HydrogenHardcoded>();
	else if (hOption == "hff2")
		_hData = std::make_unique<HydrogenFromFiles>(2);
	else if (hOption == "hff4")
		_hData = std::make_unique<HydrogenFromFiles>(4);
	else
		_hData = std::make_unique<HydrogenFromFiles>();

	// H2 data
	if (h2Option.empty())
		_h2Data = std::make_unique<H2Data>();
	else
	{
		int maxJ, maxV;
		std::istringstream() >> maxJ >> maxV;
		if (maxJ < 0 || maxV < 0)
			Error::runtime("moleculeChoice is not of a correct format. It "
				       "should be "
				       "\"maxJ maxV\"");
		_h2Data = std::make_unique<H2Data>(maxJ, maxV);
	}
}
