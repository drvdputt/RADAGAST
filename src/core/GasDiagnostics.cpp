#include "GasDiagnostics.hpp"
#include "TemplatedUtils.hpp"
#include <iostream>

namespace RADAGAST
{
    void GasDiagnostics::setHeating(const std::string& key, double value) { _heatingm[key] = value; }

    void GasDiagnostics::setCooling(const std::string& key, double value) { _coolingm[key] = value; }

    void GasDiagnostics::setUserValue(const std::string& key, double value)
    {
        auto valueIt = _userValuem.find(key);
        if (valueIt != _userValuem.end())
        {
            std::cerr << "Warning: userValues already contains value with this key. "
                         "Overwriting.\n";
            valueIt->second = value;
        }
        else
            _userValuem[key] = value;
    }
}
