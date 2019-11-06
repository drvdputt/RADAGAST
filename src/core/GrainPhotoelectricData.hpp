#ifndef CORE_GRAINPHOTOELECTRICDATA_HPP
#define CORE_GRAINPHOTOELECTRICDATA_HPP

#include "GrainPhotoelectricCalculator.hpp"

class GrainPhotoelectricData
{
public:
	/** Create new grain photoelectric data object. Takes one argument: true if
	    carbonaceous, false if silicate. The workfunction is set, and a
	    GrainPhotoelectricCalculator will be provided, which will implement the recipe of
	    Weingartner and Draine (2001). When the calculator is created, it will cache bunch
	    of things based on the given list of sizes, see its documentation. 

	    If a different photoelectric data + calculation recipe would be implemented, then
	    the abstraction could be made as follows:

	    - Make both GrainPhotoelectricData and GrainPhotoelectricCalculator abstract.

	    - makeCalculator becomes strictly virtual

	    - A concrete GrainPhotoelectricData instance will implement makeCalculator, and in
              that way polymorphically construct the right GrainPhotoelectricCalculator instance
              (Abstract Factory pattern).

	    But for now, everything is kept concrete until we implement something else than the
	    WD01 carbonaceous / silicate recipes.*/
	GrainPhotoelectricData(bool carOrSil);
	GrainPhotoelectricCalculator makeCalculator(const Array& av);
private:
	bool _carOrSil;
	double _workFunction;
};

#endif // CORE_GRAINPHOTOELECTRICDATA_HPP