#ifndef _SRC_ERROR_H_
#define _SRC_ERROR_H_

#include <iostream>

namespace Error
{

static inline void runtime(std::string message)
{
	std::cerr << "Runtime error: " << message << std::endl;
	abort();
}

template<typename T>
void rangeCheck(std::string variable, T value, T min, T max)
{
	if (value <= min || value >= max)
	{
		std::cerr << "Range error: " << variable << " = " << value << ". Should be between "
		          << min << " and " << max << std::endl;
		abort();
	}
}

template<typename T>
void reportOverridden(std::string name, T original, T replacement,
                                    std::string reason)
{
	std::cout << std::scientific << name << "has been overridden from " << original << " to "
	          << replacement << " because " << reason << std::endl;
}
} /* namespace Error */
#endif /* _SRC_ERROR_H_ */
