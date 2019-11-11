#include "GrainPhotoelectricData.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "WeingartnerDraine2001.hpp"

GrainPhotoelectricData::GrainPhotoelectricData(bool carOrSil)
                : _carOrSil{carOrSil}, _workFunction{WD01::workFunction(carOrSil)}
{
}

std::unique_ptr<GrainPhotoelectricCalculator>
GrainPhotoelectricData::makeCalculator(const Array& av) const
{
	return std::make_unique<GrainPhotoelectricCalculator>(av, _workFunction, _carOrSil);
}
