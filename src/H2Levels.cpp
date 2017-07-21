#include "H2Levels.h"
#include "H2FromFiles.h"

H2Levels::H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv) :
		NLevel(hff, frequencyv), _hff(hff)
{
}

H2Levels::~H2Levels() = default;

double H2Levels::dissociationRate(const NLevel::Solution& s) const
{
	// TODO
	return 0;
}

double H2Levels::dissociationHeating(const NLevel::Solution& s) const
{
	// TODO
	return 0;
}
