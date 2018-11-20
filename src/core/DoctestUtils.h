#ifndef GASMODULE_GIT_SRC_DOCTESTUTILS_H_
#define GASMODULE_GIT_SRC_DOCTESTUTILS_H_

#include "TemplatedUtils.h"

namespace DoctestUtils
{

inline void checkRange(std::string quantityName, double value, double min, double max)
{
	bool inRange = TemplatedUtils::inRange(value, min, max);
	CHECK_MESSAGE(inRange, quantityName << " = " << value << ". Should be between " << min
	                                    << " and " << max);
}

inline void checkTolerance(std::string quantityName, double value, double reference,
                           double precision)
{
	bool withinTolerance = TemplatedUtils::equalWithinTolerance<double>(value, reference,
	                                                                    precision);
	CHECK_MESSAGE(withinTolerance,
	              quantityName << " " << value << " is not within " << precision
	                           << " precision of reference value " << reference);
}

} // namespace DoctestUtils

#endif /* GASMODULE_GIT_SRC_DOCTESTUTILS_H_ */
