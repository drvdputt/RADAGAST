#ifndef GASMODULE_GIT_SRC_GRAININFO_H_
#define GASMODULE_GIT_SRC_GRAININFO_H_

#include <vector>

namespace GasModule
{
/** Class that a client code should use to pass the grain properties in a cell. The AbsQvv
    members should contain the absorption efficiency for each grain size (first index) and each
    frequency/wavelength. */
class GrainInfo
{
public:
	/* Well this is over engineered ... */
	enum class AvailableContents
	{
		// All members empty
		EMPTY,
		// Carbon members filled in, silicate data empty
		CARBON,
		// Vice versa
		SILICA,
		// All members filled in
		BOTH
	};

	/** Creates and empty GrainInfo. */
	GrainInfo();

	/** Constructor for carbonaceous or silicate only. */
	enum class GrainType
	{
		CARBON,
		SILICA
	};
	/** Converts the grain type argument of this constructor to a status. */
	static AvailableContents SingleTypeAvailable(GrainType t)
	{
		return t == GrainType::CARBON ? AvailableContents::CARBON
		                              : AvailableContents::SILICA;
	}

	GrainInfo(GrainType t, const std::vector<double>& grainSizev,
	          const std::vector<std::vector<double>>& absQvv);

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInfo(const std::vector<double>& carbonGrainSizev,
	          const std::vector<std::vector<double>>& carbonAbsQvv,
	          const std::vector<double>& silicaGrainSizev,
	          const std::vector<std::vector<double>>& silicaAbsQvv);

private:
	const AvailableContents _contents{AvailableContents::EMPTY};

	typedef struct GrainData
	{
		GrainData() = default;
		GrainData(const std::vector<double>& grainSizev,
		          const std::vector<std::vector<double>> absQvv)
		                : _grainSizev{grainSizev}, _absQvv{absQvv}
		{
		}
		const std::vector<double> _grainSizev;
		const std::vector<std::vector<double>> _absQvv;
	} GrainData;

	// Carbonaceous (graphite and PAH) grains
	const GrainData _carbon;
	// Silicate grains
	const GrainData _silica;
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININFO_H_ */
