#include "GrainPhotoelectricData.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "WeingartnerDraine2001.hpp"

namespace GasModule
{
    GrainPhotoelectricData::GrainPhotoelectricData(bool carOrSil)
        : _carOrSil{carOrSil}, _workFunction{WD01::workFunction(carOrSil)}
    {}

    std::unique_ptr<GrainPhotoelectricCalculator>
    GrainPhotoelectricData::makeCalculator(const Array* av, const std::vector<Array>* qAbsvv,
                                           const Spectrum* meanIntensity) const
    {
        return std::make_unique<GrainPhotoelectricCalculator>(av, qAbsvv, _workFunction, _carOrSil, meanIntensity);
    }
}
