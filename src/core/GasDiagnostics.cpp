#include "GasDiagnostics.hpp"
#include "TemplatedUtils.hpp"
#include <iostream>

namespace RADAGAST
{
    void GasDiagnostics::setReactionInfo(const std::vector<std::string>& reactionNamev,
                                         const std::vector<double>& reactionRatev)
    {
        for (size_t i = 0; i != reactionNamev.size(); ++i)
            _reactionInfov["chem:" + reactionNamev[i]] = reactionRatev[i];
    }

    void GasDiagnostics::setHeating(const std::string& name, double value) { _heatingm["heat:" + name] = value; }

    void GasDiagnostics::setCooling(const std::string& name, double value) { _coolingm["cool:" + name] = value; }

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
