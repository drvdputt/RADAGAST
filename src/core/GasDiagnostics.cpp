#include "GasDiagnostics.h"
#include "TemplatedUtils.h"

#include <iostream>

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
