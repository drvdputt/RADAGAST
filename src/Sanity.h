#ifndef _SRC_SANITY_H_
#define _SRC_SANITY_H_

#include "flags.h"
#include <ios>
#include <iostream>
#include <string>

namespace Sanity
{
void reportOverridden(std::string name, double original, double replacement, std::string reason)
{
	std::cout << std::scientific << name << "has been overridden from " << original << " to " << replacement << " because "
			<< reason << std::endl;
}
} /* namespace Sanity */
#endif /* _SRC_SANITY_H_ */
