#ifndef _TWOLEVEL_H_
#define _TWOLEVEL_H_

#include "NLevel.h"

class TwoLevel : public NLevel
{
public:
	/* Creates an object that represents a two-level (component of a) medium. The level
	   population equilibrium will be calculated using the frequency grid supplied as argument
	   of the constructor. */
	TwoLevel();
};

#endif /* _TWOLEVEL_H_ */
