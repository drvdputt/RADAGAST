#ifndef _SRC_ERROR_H_
#define _SRC_ERROR_H_

#include "global.h"
#include <ios>
#include <iostream>
#include <string>

namespace Error
{

static inline void runtime(std::string message)
{
	std::cerr << "Runtime error: " << message << std::endl;
	abort();
}

static inline void rangeCheck(std::string variable, double value, double min, double max)
{
	if (value < min || value > max)
	{
		std::cerr << "Range error: " << variable << " = " << value << ". Should be between "
		          << min << " and " << max << std::endl;
		abort();
	}
}

static inline void reportOverridden(std::string name, double original, double replacement,
                                    std::string reason)
{
	std::cout << std::scientific << name << "has been overridden from " << original << " to "
	          << replacement << " because " << reason << std::endl;
}
} /* namespace Error */
#endif /* _SRC_ERROR_H_ */
