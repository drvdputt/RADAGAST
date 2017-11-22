#ifndef GASMODULE_GIT_SRC_HYDROGENDATAPROVIDER_H_
#define GASMODULE_GIT_SRC_HYDROGENDATAPROVIDER_H_

#include "LevelDataProvider.h"

#include <array>

/** This class exists to underscore the main difference between a regular LevelDataProvider and one
    suited for hydrogen. Using a LevelDataProvider of this subclass the two-photon continuum between
    1s and 2s can be implemented, as it provides a function to retrieve the right indices. The rest
    of the details about the levels remains hidden. */
class HydrogenDataProvider : public LevelDataProvider
{
protected:
	HydrogenDataProvider() = default;
	virtual ~HydrogenDataProvider() = default;
public:
	/** Returns the upper and lower index (in that order) of the two photon transition as an
	    array of size 2. This is a pure virtual function, as the indices returned depend on the
	    indexing scheme used, which depends on the specific implementation of how the hydrogen
	    data is loaded. Therefore, subclasses which use different data will have different
	    implementations for this. */
	virtual std::array<size_t, 2> twoPhotonIndices() const = 0;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENDATAPROVIDER_H_ */
