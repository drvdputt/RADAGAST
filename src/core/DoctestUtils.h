#ifndef GASMODULE_GIT_SRC_DOCTESTUTILS_H_
#define GASMODULE_GIT_SRC_DOCTESTUTILS_H_

#include "doctest.h"
#include "TemplatedUtils.h"

#include <sstream>

/** Here are some wrappers around some of the doctest assertion macros. Frequently used checks
    are defined here, with some standardized, readable failure messages. */
namespace DoctestUtils
{

inline void checkRange(std::string quantityName, double value, double min, double max,
                       bool warn = false)
{
	bool inRange = TemplatedUtils::inRange(value, min, max);
	if (!inRange)
	{
		std::stringstream ss;
		ss << quantityName << " = " << value << ". Should be between " << min << " and "
		   << max;
		std::string message = ss.str();
		if (warn)
			WARN_MESSAGE(inRange, message);
		else
			CHECK_MESSAGE(inRange, message);
	}
}

inline void checkTolerance(std::string quantityName, double value, double reference,
                           double precision, bool warn = false)
{
	bool withinTolerance = TemplatedUtils::equalWithinTolerance<double>(value, reference,
	                                                                    precision);
	if (!withinTolerance)
	{
		std::stringstream ss;
		ss << quantityName << " " << value << " is not within " << precision
		   << " precision of reference value " << reference;
		std::string message = ss.str();
		if (warn)
			WARN_MESSAGE(withinTolerance, message);
		else
			CHECK_MESSAGE(withinTolerance, message);
	}
}

} // namespace DoctestUtils

#endif /* GASMODULE_GIT_SRC_DOCTESTUTILS_H_ */
