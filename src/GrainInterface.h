#ifndef GASMODULE_GIT_SRC_GRAININTERFACE_H_
#define GASMODULE_GIT_SRC_GRAININTERFACE_H_

#include <valarray>
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
    are listed in the nested 'Population' class. Another class, closely related to this one but not
    to be included by the client code, is GrainTypeProperties. It contains a bunch of functions that
    return parameters pertaining to the different choices of GrainType listed above. Those functions
    are in a separate file because they have no use in the public interface. */
class GrainInterface
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
		Population(GrainType type, const std::valarray<double>& sizev,
		           const std::valarray<double>& densityv,
		           const std::valarray<double>& temperaturev,
		           const std::vector<std::valarray<double>>& qAbsvv);

		const GrainType _type;
		const std::valarray<double> _sizev, _densityv, _temperaturev;
		const std::vector<std::valarray<double>> _qAbsvv;
	};

	/** Creates and empty GrainInfo. */
	GrainInterface();

	/** Constructor for a mix of carbonaceous and silicate. */
	GrainInterface(const std::valarray<double>& carbonaceousGrainSizev,
	               const std::valarray<double>& carbonaceousDensityv,
	               const std::vector<std::valarray<double>>& carbonaceousAbsQvv,
	               const std::valarray<double>& carTv,
	               const std::valarray<double>& silicateGrainSizev,
	               const std::valarray<double>& silicateDensityv,
	               const std::vector<std::valarray<double>>& silicateAbsQvv,
	               const std::valarray<double>& silTv);

	/** Constructor which takes a vector of predefined populations. */
	GrainInterface(const std::vector<Population>& populationv);

	std::vector<Population> populationv() const { return _populationv; }

private:
	const std::vector<Population> _populationv{};
};
} /* namespace GasModule */

#endif /* GASMODULE_GIT_SRC_GRAININTERFACE_H_ */
