#ifndef CORE_HDATA_HPP
#define CORE_HDATA_HPP

#include "LevelCoefficients.hpp"
#include "RecombinationRate.hpp"

#include <array>

/** Abstract class which specifies how to provide data for the hydrogen model. It inherits from
    LevelCoefficients, so that level coefficients can be delivered in a consistent way. The
    recombination coefficients are concretely implemented here. */
class HData : public LevelCoefficients
{
public:
	HData();
	virtual ~HData() = default;

	/** The maximum quantum number n for which data is available */
	virtual int nMax() const = 0;

	/** Get the level index for a given pair of quantum numbers n, l. If the data is not
	    resolved on l, the l argument is ignored. */
	virtual size_t index(int n, int l) const = 0;

	/** Returns the upper and lower index (in that order) of the two photon transition as an
	    array of size 2. This is a pure virtual function, as the indices returned depend on the
	    indexing scheme used, which depends on the specific implementation of how the hydrogen
	    data is loaded. Therefore, subclasses which use different data will have different
	    implementations for this. */
	virtual std::array<size_t, 2> twoPhotonIndices() const = 0;

	/** Get the recombination rate to each level. Multiply with ne * np to get cm-3 s-1.
	    [cm3 s-1] */
	EVector recombinationRatev(double T) const;

private:
	std::unique_ptr<const RecombinationRate> _rr;
};

#endif // CORE_HDATA_HPP
