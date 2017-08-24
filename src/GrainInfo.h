#ifndef GASMODULE_GIT_SRC_GRAININFO_H_
#define GASMODULE_GIT_SRC_GRAININFO_H_

#include "Array.h"

#include <vector>

namespace GasModule
{
/** List of possible grain types */
enum class GrainType
{
	CAR,
	SIL
};

/** Class that a client code should use to pass the grain properties in a cell. This class contains
    a list of grain models, for each of which a set of properties must be given. These properties,
    are listed in the nested 'Population' class. */
class GrainInfo
{
public:
	/** The properties that need to be given per grain model that needs to be included. The
	    grain type can be one of the types listed in the enum above. Then, an array of sizes,
	    number densities, and one of temperatures needs to be specified (consider even using
	    grain temperature distributions (for each grain size!) in the future). The temperatures
	    are important for the H2 formation rate on the surfaces of the grains. The absorption
	    efficiency also needs to be given (need to review its exact definition). It represents
	    the number of photons absorbed from an ambient radiation field at a certain wavelength.
	    It needs to be given for each grain size, and for each point of the frequency grid. */
	class Population
	{
	public:
		Population(GrainType type, const Array& sizev, const Array& densityv,
		           const Array& temperaturev, const std::vector<Array>& qAbsvv);

		const GrainType _type;
		const Array _sizev, _densityv, _temperaturev;
		const std::vector<Array> _qAbsvv;
	};

	/** Creates and empty GrainInfo. */
	GrainInfo();

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInfo(const Array& carbonaceousGrainSizev, const Array& carbonaceousDensityv,
	          const std::vector<Array>& carbonaceousAbsQvv, const Array& carTv,
	          const Array& silicateGrainSizev, const Array& silicateDensityv,
	          const std::vector<Array>& silicateAbsQvv, const Array& silTv);

	std::vector<Population> populationv() const { return _populationv; }

private:
	const std::vector<Population> _populationv{};
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININFO_H_ */
