#ifndef _NLEVEL_H_
#define _NLEVEL_H_

#include "Array.h"

class NLevel
{
public:
	NLevel(const Array& frequencyv);

	void solveBalance(double atomDensity, double electronDensity, double protonDensity, double T, const Array& specificIntensityv);

	Array emissivityv() const;

	Array opacityv() const;

	Array scatteringOpacityv() const;
	Array absorptionOpacityv() const;

private:

};

#endif /* _NLEVEL_H_ */
