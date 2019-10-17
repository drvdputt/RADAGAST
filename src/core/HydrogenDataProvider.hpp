#ifndef CORE_HYDROGENDATAPROVIDER_H_
#define CORE_HYDROGENDATAPROVIDER_H_

#include "Constants.hpp"
#include "LevelDataProvider.hpp"
#include "NLevel.hpp"
#include "RecombinationRate.hpp"

#include <array>

/** Abstract class which specifies how to provide data for the hydrogen model. It inherits from
    LevelCoefficients, so that level coefficients can be delivered in a consistent way. */
class HData : public LevelCoefficients
{
public:
	HData()
	                : LevelCoefficients(Constant::HMASS_CGS),
	                  _rr{std::make_unique<HydrogenADF48>()} {};
	virtual ~HData();

	virtual int nMax() const = 0;
	virtual size_t index(int n, int l) const = 0;

	/** Returns the upper and lower index (in that order) of the two photon transition as an
	    array of size 2. This is a pure virtual function, as the indices returned depend on the
	    indexing scheme used, which depends on the specific implementation of how the hydrogen
	    data is loaded. Therefore, subclasses which use different data will have different
	    implementations for this. */
	virtual std::array<size_t, 2> twoPhotonIndices() const = 0;

private:
	std::unique_ptr<const RecombinationRate> _rr;
};

#endif /* CORE_HYDROGENDATAPROVIDER_H_ */
