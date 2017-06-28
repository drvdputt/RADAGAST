#ifndef _SRC_HYDROGENLEVELS_H_
#define _SRC_HYDROGENLEVELS_H_

#include "NLevel.h"

class HydrogenLevels : public NLevel
{
public:
	HydrogenLevels(bool hardcoded);

	Array emissivityv(const Solution& s) const override;

private:
	LevelDataProvider* chooseDataProvider(bool hardcoded) const;

protected:
	Array boundBoundContinuum(const Solution& s) const;
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
