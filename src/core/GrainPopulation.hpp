#ifndef CORE_GRAINPOPULATION_HPP
#define CORE_GRAINPOPULATION_HPP

#include "Array.hpp"

#include <vector>

/** List of grain types which have built-in values for H2 formation and the photoelectric
    effect. For anything else, the label 'OTHER' can be used, but then no contribution to the H2
    formation or photoelectric heating will be made by this grain population. */
enum class GrainTypeLabel
{
	CAR,
	SIL,
	OTHER
};

class GrainH2Formation;
class GrainPhotoelectricData;

class GrainPopulation
{
public:
	/** The grain type can be one of the types listed in the enum above. Then, an array of
            sizes, number densities, and one of temperatures needs to be specified (consider
            even using grain temperature distributions (for each grain size!) in the future).
            The temperatures are important for the H2 formation rate on the surfaces of the
            grains. The absorption efficiency also needs to be given. It needs to be given for
            each grain size, and for each point of the frequency grid of the input radiation
            field. */
	GrainPopulation(GrainTypeLabel type, const Array& sizev, const Array& densityv,
	                const Array& temperaturev, const std::vector<Array>& qAbsvv);

	/** Undelete the move constructor (needed to be able to put these objects into a vector
	    std::vector) */
	GrainPopulation(GrainPopulation&&);
	~GrainPopulation();

	/** @name Trivial getters. */
	/**@{*/
	std::vector<Array> qAbsvv() const { return _qAbsvv; }
	std::valarray<double> temperaturev() const { return _temperaturev; }
	std::valarray<double> densityv() const { return _densityv; }
	std::valarray<double> sizev() const { return _sizev; }
	/**@}*/

	size_t numSizes() const { return _sizev.size(); }

	/** If not a nullptr, then the object can be used to calculate the h2 formation rate. */
	const GrainH2Formation* h2formationData() { return _h2formation.get(); }

	/** If not a nullptr, then the object can be used to calculate the photoelectric effect
	    (charge distribution and heating rate). */
	const GrainPhotoelectricData* photoelectricData() { return _photoelectricData.get(); }

	Array qAbsv(int m) const { return _qAbsvv[m]; }
	double temperature(int m) const { return _temperaturev[m]; }
	double density(int m) const { return _densityv[m]; }
	double size(int m) const { return _sizev[m]; }

	void test() const;

private:
	Array _sizev;
	Array _densityv;
	Array _temperaturev;
	std::vector<Array> _qAbsvv;
	std::unique_ptr<GrainH2Formation> _h2formation{nullptr};
	std::unique_ptr<GrainPhotoelectricData> _photoelectricData{nullptr};
};

#endif // CORE_GRAINPOPULATION_HPP
