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
		CAR,
		SIL
	};
	GrainInfo(GrainType t, const Array& grainSizev, const Array& grainDensityv,
	          const std::vector<Array>& absQvv, const Array& Tv);

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInfo(const Array& carbonaceousGrainSizev, const Array& carbonaceousDensityv,
	          const std::vector<Array>& carbonaceousAbsQvv, const Array& carTv,
	          const Array& silicateGrainSizev, const Array& silicateDensityv,
	          const std::vector<Array>& silicateAbsQvv, const Array& silTv);
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
		                      const std::vector<Array> absQvv, const Array& Tv,
		                      const GrainInfo::SurfaceInteractionParameters& par)
		                : _grainSizev{grainSizev},
		                  _grainDensityv{grainDensityv}, _absQvv{absQvv}, _tv{Tv}, _par{par}
		{
		}
		const Array _grainSizev{};
		const Array _grainDensityv{};
		const std::vector<Array> _absQvv{};
		const Array _tv{};
		const GrainInfo::SurfaceInteractionParameters _par;
	} GrainSizeDistribution;

	// Carbonaceous (graphite and PAH) grains
	const GrainSizeDistribution _carSizeDist{};
	// Silicate grains
	const GrainSizeDistribution _silSizeDist{};
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININFO_H_ */
