#ifndef _SRC_ERROR_H_
#define _SRC_ERROR_H_

#include <iostream>

namespace Error
{

/** Prints a message to stderr and aborts. */
static inline void runtime(std::string message)
{
	std::cerr << "Runtime error: " << message << std::endl;
	abort();
}

template <typename T> void equalCheck(std::string variable_names, T value1, T value2)
{
	if (value1 != value2)
	{
		std::cerr << "Equality error: " << variable_names
		          << " should be equal. Their values are " << value1 << " and " << value2
		          << std::endl;
		abort();
	}
}

/** Prints a message to stderr and aborts if the given variable does not lie between min and max. */
template <typename T> void rangeCheck(std::string variable, T value, T min, T max)
{
	if (value <= min || value >= max)
	{
		std::cerr << "Range error: " << variable << " = " << value << ". Should be between "
		          << min << " and " << max << std::endl;
		abort();
	}
}

/** Can probably be removed */
template <typename T>
void reportOverridden(std::string name, T original, T replacement, std::string reason)
{
	std::cout << std::scientific << name << "has been overridden from " << original << " to "
	          << replacement << " because " << reason << std::endl;
}
} /* namespace Error */
#endif /* _SRC_ERROR_H_ */
