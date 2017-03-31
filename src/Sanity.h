#ifndef _SRC_SANITY_H_
#define _SRC_SANITY_H_

#include <ios>
#include <iostream>
#include <string>
#include "global.h"

namespace Sanity
{
void reportOverridden(std::string name, double original, double replacement, std::string reason)
{
	std::cout << std::scientific << name << "has been overridden from " << original << " to "
	          << replacement << " because " << reason << std::endl;
}
} /* namespace Sanity */
#endif /* _SRC_SANITY_H_ */
