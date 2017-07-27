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
		CARBON,
		SILICA
	};
	GrainInfo(GrainType t, const Array& grainSizev,
	          const std::vector<Array>& absQvv);

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInfo(const Array& carbonGrainSizev,
	          const std::vector<Array>& carbonAbsQvv,
	          const Array& silicaGrainSizev,
	          const std::vector<Array>& silicaAbsQvv);

	bool hasCarbon() const {return _hasCarbon;}
	bool hasSilica() const {return _hasSilica;}


private:
	bool _hasCarbon{false};
	bool _hasSilica{false};

public:
	typedef struct GrainSizeDistribution
	{
		GrainSizeDistribution() = default;
		GrainSizeDistribution(const Array& grainSizev,
		          const std::vector<Array> absQvv)
		                : _grainSizev{grainSizev}, _absQvv{absQvv}
		{
		}
		const Array _grainSizev;
		const Array _grainDensityv;
		const std::vector<Array> _absQvv;
	} GrainData;

	// Carbonaceous (graphite and PAH) grains
	const GrainSizeDistribution _carbonSizeDist;
	// Silicate grains
	const GrainSizeDistribution _silicaSizeDist;
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININFO_H_ */
