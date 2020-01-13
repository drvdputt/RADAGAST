#ifndef CORE_ERROR_HPP
#define CORE_ERROR_HPP

#include <iostream>

namespace Error
{
    static inline void runtime(std::string message);
    template<typename T> void equalCheck(std::string variable_names, T value1, T value2);
    template<typename T> void rangeCheck(std::string variable, T value, T min, T max);
}

#include "TemplatedUtils.hpp"
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
    template<typename T> void equalCheck(std::string variable_names, T value1, T value2)
    {
        if (value1 != value2)
        {
            std::cerr << "Equality error: " << variable_names << " should be equal. Their values are " << value1
                      << " and " << value2 << std::endl;
            abort();
        }
    }

    /** Prints a message to stderr and aborts if the given variable does not lie between min and
    max. */
    template<typename T> void rangeCheck(std::string variable, T value, T min, T max)
    {
        if (!TemplatedUtils::inRange(value, min, max))
        {
            std::cerr << "Range error: " << variable << " = " << value << ". Should be between " << min << " and "
                      << max << std::endl;
            abort();
        }
    }

} /* namespace Error */
#endif  // CORE_ERROR_HPP
