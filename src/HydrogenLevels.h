#ifndef _SRC_HYDROGENLEVELS_H_
#define _SRC_HYDROGENLEVELS_H_

#include "HydrogenDataProvider.h"
#include "NLevel.h"

class HydrogenLevels : public NLevel
{
public:
	// TODO: consider move semantics here?
	/* Work with a shared pointer here, to make construction using a temporary object
	   possible. */
	HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp);

	Array emissivityv(const Solution& s) const override;

private:
	Array twoPhotonEmissivityv(const Solution& s) const;

	std::shared_ptr<const HydrogenDataProvider> _hdp;
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
