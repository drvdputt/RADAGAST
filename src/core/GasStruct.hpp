#ifndef CORE_GASSTRUCT_H_
#define CORE_GASSTRUCT_H_

#include "EigenAliases.hpp"
#include "SpeciesIndex.hpp"

/** A struct containing a bunch of parameters about the gas that are frequently passed around,
    such as temperature, the density vector, and the specific intensity of the ambient radiation
    field. Can probably be expanded without too much overhead.

    This struct is not to be confused with GasState. Objects of the GasState class give the
    client code a way to store the result of the calculation.

    This struct has a distinctly different purpose, which is passing around gas properties to
    the different functions of the code. If we decide that an extra property is needed somewhere
    in a subclass' specialization of some member function, we can simple expand this object, instead of having to
    change the function signature for the parent class and all subclasses.

    I personally prefer that this object should always be passed around as a const reference.
    Letting functions modify this blob of parameters seems like a bad idea for clarity to me.
    Ideally functions only take as an argument what they need to read, and then return a new
    version of whatever they've modified.

    So, in other words, just to be clear: if you want to modify an environment parameter, just
    do something like

    environment.p = function_giving_new_p(const environment&);

    Of course, I have no way to enforce this. Making the members constant would make even the
    above impossible. Just consider it the style of the code to avoid passing by non-const
    reference as much as possible. */
typedef struct GasStruct
{
	/** Default constructor, which just provides safe parameters for a vacuum. T - 0 can
	    create nans, so we pick a random constant T here. */
	GasStruct() : _T{500}, _speciesNv{EVector::Zero(SpeciesIndex::size())} {}

	/** Package the temperature and the species vector into one objects, as these are
	    typically needed together. */
	GasStruct(double T, const EVector& speciesNv) : _T{T}, _speciesNv{speciesNv} {}

	double _T;
	EVector _speciesNv;

	// The fraction of ortho H2 might be important too. Default it here to .25, which
	// corresponds to a ratio of 3 to 1
	double _orthoH2{.75};
} GasStruct;

#endif /* CORE_GASSTRUCT_H_ */
