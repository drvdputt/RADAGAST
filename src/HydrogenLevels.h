#ifndef GASMODULE_GIT_SRC_HYDROGENLEVELS_H_
#define GASMODULE_GIT_SRC_HYDROGENLEVELS_H_

#include "NLevel.h"

class HydrogenDataProvider;

class HydrogenLevels: public NLevel
{
public:
	/** Work with a shared pointer here, so that a pointer can be given both to the derived
	    class and the base class. */
	HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp, const Array& frequencyv);
	~HydrogenLevels();

	Array emissivityv(const Solution& s) const override;

private:
	Array twoPhotonEmissivityv(const Solution& s) const;

	std::shared_ptr<const HydrogenDataProvider> _hdp;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENLEVELS_H_ */
