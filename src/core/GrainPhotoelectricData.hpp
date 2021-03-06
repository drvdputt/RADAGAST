#ifndef CORE_GRAINPHOTOELECTRICDATA_HPP
#define CORE_GRAINPHOTOELECTRICDATA_HPP

#include "GrainPhotoelectricCalculator.hpp"

namespace RADAGAST
{
    /** This class stores data related to the photoelectric effect, which is static (i.e. does
        not depend on grain sizes, no caching). It also acts as a factory for
        GrainPhotoelectricCalculator objects. If a different photoelectric data + calculation
        recipe would be implemented, then the abstraction could be made as follows (this
        approach is probably called abstract factory):

        - Make both GrainPhotoelectricData and GrainPhotoelectricCalculator abstract

        - makeCalculator becomes strictly virtual

        - A concrete GrainPhotoelectricData instance will implement makeCalculator, and in that way
          polymorphically construct the right GrainPhotoelectricCalculator instance (Abstract
          Factory pattern).

          But for now, everything is kept concrete until we implement something else than the WD01
          carbonaceous / silicate recipes.*/
    class GrainPhotoelectricData
    {
    public:
        /** Create new grain photoelectric data object. Takes one argument: true if carbonaceous,
            false if silicate. The workfunction is set, and a GrainPhotoelectricCalculator will be
            provided, which will implement the recipe of Weingartner and Draine (2001). When the
            calculator is created, it will cache bunch of things based on the given list of sizes,
            see its documentation. */
        GrainPhotoelectricData(bool carOrSil);

        /** Create a photoelectric calculator, for a given array of grain sizes, a constant
            radiation field, and Qabs values for each size and frequency. This object will serve
            as a workspace for the photoelectric heating, allowing the storage of intermediate
            results and caching. Passing the requested constants here helps with some
            precalculations and caching mechanisms inside the GrainPhotoelectricCalculator. */
        std::unique_ptr<GrainPhotoelectricCalculator>
        makeCalculator(const Array* sizev, const std::vector<Array>* qAbsvv, const Spectrum* meanIntensity) const;

    private:
        bool _carOrSil;
        double _workFunction;
    };
}
#endif  // CORE_GRAINPHOTOELECTRICDATA_HPP
