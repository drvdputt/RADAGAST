#include "HydrogenLevels.hpp"
#include "Constants.hpp"
#include "GasStruct.hpp"
#include "HydrogenDataProvider.hpp"
#include "IonizationBalance.hpp"
#include "Options.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"

#include <vector>

using namespace std;

Array HydrogenLevels::emissivityv(const Solution& s, const Array& eFrequencyv) const
{
	return lineEmissivityv(s, eFrequencyv) + twoPhotonEmissivityv(s, eFrequencyv);
}

