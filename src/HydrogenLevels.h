#ifndef _SRC_HYDROGENLEVELS_H_
#define _SRC_HYDROGENLEVELS_H_

#include "NLevel.h"
#include "Table.h"

#include <vector>

class HydrogenLevels : public NLevel
{
public:
	HydrogenLevels(bool hardcoded);

private:
	LevelDataProvider* chooseDataProvider(bool hardcoded) const;

protected:
	Array boundBoundContinuum(const Solution& s) const override;
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
