#ifndef GASMODULE_GIT_SRC_GRAININFO_H_
#define GASMODULE_GIT_SRC_GRAININFO_H_

#include "Array.h"

#include <vector>

namespace GasModule
{
/** Class that a client code should use to pass the grain properties in a cell. The AbsQvv
    members should contain the absorption efficiency for each grain size (first index) and each
    frequency/wavelength. */
class GrainInfo
{
public:
	/** Creates and empty GrainInfo. */
	GrainInfo();

	/** Constructor for carbonaceous or silicate only. */
	enum class GrainType
	{
		CARBONACEOUS,
		SILICATE
	};
	GrainInfo(GrainType t, const Array& grainSizev, const Array& grainDensityv,
	          const std::vector<Array>& absQvv);

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInfo(const Array& carbonGrainSizev, const Array& carbonaceousDensityv,
	          const std::vector<Array>& carbonAbsQvv, const Array& silicaGrainSizev,
	          const Array& silicateDensityv, const std::vector<Array>& silicaAbsQvv);

	bool hasCarbonaceous() const { return _hasCar; }
	bool hasSilicate() const { return _hasSil; }

private:
	bool _hasCar{false};
	bool _hasSil{false};

public:
	typedef struct GrainSizeDistribution
	{
		GrainSizeDistribution() = default;
		GrainSizeDistribution(const Array& grainSizev, const Array& grainDensityv,
		                      const std::vector<Array> absQvv)
		                : _grainSizev{grainSizev},
		                  _grainDensityv{grainDensityv}, _absQvv{absQvv}
		{
		}
		const Array _grainSizev;
		const Array _grainDensityv;
		const std::vector<Array> _absQvv;
	} GrainSizeDistribution;

	// Carbonaceous (graphite and PAH) grains
	const GrainSizeDistribution _carSizeDist{};
	// Silicate grains
	const GrainSizeDistribution _silSizeDist{};
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININFO_H_ */
