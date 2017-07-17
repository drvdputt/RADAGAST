#include "H2Levels.h"
#include "H2FromFiles.h"

H2Levels::H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv) :
		NLevel(hff, frequencyv), _hff(hff)
{
}

H2Levels::~H2Levels() = default;
