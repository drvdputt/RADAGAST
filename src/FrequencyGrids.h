#include "Array.h"

/** Collection of the different frequency grids used in the code */
class FrequencyGrids
{
public:
	FrequencyGrid(const Array& iFrequencyv, const Array& eFrequencyv,
	              const Array& oFrequencyv)
	                : _iFrequencyv(iFrequencyv), _eFrequencyv(eFrequencyv),
	                  _oFrequencyv(oFrequencyv)
	{
	}

	const Array _iFrequencyv, _eFrequencyv, _oFrequencyv;
};
