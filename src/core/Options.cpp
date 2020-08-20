#include "Options.hpp"

namespace RADAGAST
{

    std::string Options::getenvWithDefault(const std::string& variableName, const std::string& defaultValue)
    {
        const char* value = std::getenv(variableName.c_str());
        return value ? value : defaultValue;
    }

}
