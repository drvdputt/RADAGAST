#ifndef GASMODULE_GIT_SRC_ERROR_H_
#define GASMODULE_GIT_SRC_ERROR_H_

#include <iostream>

namespace Error
{
	static inline void runtime(std::string message);
	template <typename T> void equalCheck(std::string variable_names, T value1, T value2);
	template <typename T> void rangeCheck(std::string variable, T value, T min, T max);
	static inline void fuzzyCheck(std::string variable, double value, double reference,
	                       double precision);
}

#include "TemplatedUtils.h"

namespace Error
{
/** Prints a message to stderr and aborts. */
static inline void runtime(std::string message)
{
	std::cerr << "Runtime error: " << message << std::endl;
	abort();
}

/** Checks if the two given values are not equal, and prints a message containing the given name
    if this is the case. */
template <typename T> void equalCheck(std::string variable_names, T value1, T value2)
{
	if (value1 != value2)
	{
		std::cerr << "Equality error: " << variable_names
		          << " should be equal. Their values are " << value1 << " and "
		          << value2 << std::endl;
		abort();
	}
}

/** Prints a message to stderr and aborts if the given variable does not lie between min and
    max. */
template <typename T> void rangeCheck(std::string variable, T value, T min, T max)
{
	if (!TemplatedUtils::inRange(value, min, max))
	{
		std::cerr << "Range error: " << variable << " = " << value
		          << ". Should be between " << min << " and " << max << std::endl;
		abort();
	}
}

/** Throws an error if @c value and @c reference are equal up to a factor @c precision */

static inline void fuzzyCheck(std::string variable, double value, double reference, double precision)
{
	if (!TemplatedUtils::equalWithinTolerance<double>(value, reference, precision))
	{
		std::cerr << variable << " " << value << " is not within" << precision
		          << " precision of reference value " << reference << std::endl;
		abort();
	}
}

} /* namespace Error */
#endif /* GASMODULE_GIT_SRC_ERROR_H_ */
