#ifndef CORE_GRAININTERFACE_HPP
#define CORE_GRAININTERFACE_HPP

#include "GrainPopulation.hpp"

#include <functional>
#include <memory>
#include <valarray>
#include <vector>

namespace GasModule
{

// NEEDED FOR H2 FORMATION //
/** Parameters for the formation of H2 on the surface on the grain. See 2013-Röllig et al. table
    C.1. */
typedef struct SfcInteractionPar
{
	SfcInteractionPar() = default;
	SfcInteractionPar(double EH2, double Es, double EHp, double EHc, double aSqrt,
	                  double nuH2, double nuHc, double F);
	/* This boolean is set to false when the default constructor is called, signifiying that
	   the created object does not contain any useful information (meaning that the
	   graintype for which it was contructed is not supported, and that the graintype given
	   in the corresponding static function should be skipped for the H2 formation rate
	   calculation. */
	bool _valid{false};
	const double _eH2{0}, _es{0}, _eHp{0}, _eHc{0}, _aSqrt{0}, _nuH2{0}, _nuHc{0}, _f{0};
} SfcInteractionPar;

/** List of grain types which have built-in values for these properties. For anything else, a
    separate graintype subclass will have to be written */
enum class GrainTypeLabel
{
	CAR,
	SIL,
};

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
	std::unique_ptr<std::vector<GrainPopulation>> _populationv;
};
} /* namespace GasModule */

#endif // CORE_GRAININTERFACE_HPP
