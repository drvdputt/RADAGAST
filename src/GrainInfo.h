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

	// See 2013-Röllig et al. table C.1.
	typedef struct SurfaceInteractionParameters
	{
		SurfaceInteractionParameters() = default;
		SurfaceInteractionParameters(double EH2, double Es, double EHp, double EHc,
		                             double aSqrt, double nuH2, double nuHc, double F)
		                : _eH2{EH2}, _es{Es}, _eHp{EHp}, _eHc{EHc}, _aSqrt{aSqrt},
		                  _nuH2{nuH2}, _nuHc{nuHc}, _f{F}
		{
		}
		const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0},
		                _f{0};
	} SurfaceInteractionParameters;

	static SurfaceInteractionParameters sfcInteractionPar(GrainType t)
	{
		if (t == GrainType::CAR)
			return {520, 260, 800, 30000, 14, 3e12, 1.3e13, 1e-10};
		else /* if  (t == GrainType::SILICATE) */
			return {320, 110, 450, 30000, 14.4, 3e12, 1.3e13, 1e-10};
	}

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
