#ifndef CORE_GRAININTERFACE_HPP
#define CORE_GRAININTERFACE_HPP

#include "GrainPopulation.hpp"

#include <functional>
#include <memory>
#include <valarray>
#include <vector>

namespace GasModule
{
/** Class that a client code should use to pass the grain properties in a cell. This class
    contains a list of grain models, for each of which a set of properties must be given. These
    properties, are listed in the nested 'Population' class. Another class, closely related to
    this one but not to be included by the client code, is GrainTypeProperties. It contains a
    bunch of functions that return parameters pertaining to the different choices of GrainType
    listed above. Those functions are in a separate file because they have no use in the public
    interface. */
class GrainInterface
{
public:
	/** Constructor which takes a vector of predefined populations. Please pass the vector
	    using a unique pointer and std::move. Ownership over the vector of populations will
	    then be transferred to this class, preventing unnecessary copying (I'm not sure what
	    the performance impact of copying would be in the future). */
	GrainInterface(std::unique_ptr<std::vector<GrainPopulation>> populationvToMove);

	/** Default constructor, equivalent to no grains at all. */
	GrainInterface();
	~GrainInterface();

	/** Undelete the move constructor, return this object from functions etc. This is
	    necessary because the presence of a unique_ptr. Very annoying if you don't know
	    about this behaviour, especially since the compiler errors can be quite cryptic. */
	GrainInterface(GrainInterface&&);

	/** The number of grain populations. */
	size_t numPopulations() const;

	/** Get a pointer to the @i 'th population. Throws an error when out of range. */
	const GrainPopulation* population(size_t i) const;

	/** Get a pointer to the whole population vector. */
	std::vector<GrainPopulation>* populationv() const;

	/** A quick test to see if all population objects in the vector have reasonable
	    contents. */
	void test() const;

private:
	std::unique_ptr<std::vector<GrainPopulation>> _populationv{nullptr};
};
} /* namespace GasModule */

#endif // CORE_GRAININTERFACE_HPP
