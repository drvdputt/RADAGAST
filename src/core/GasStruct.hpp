#ifndef CORE_GASSTRUCT_HPP
#define CORE_GASSTRUCT_HPP

#include "EigenAliases.hpp"
#include "SpeciesIndex.hpp"

/** A struct containing some things which are often needed together. This makes it possible to
    make the function signature uniform over some subclasses. This is mostly needed for the
    level calculations, where the densities, temperature and maybe ortho-to-para ratio are used
    to calculate collision coefficients. TODO: Rework and at least rename this to somthing less
    ambiguous. */
typedef struct GasStruct
{
	GasStruct(double t, const SpeciesVector& sv) : _t{t}, _sv{sv} {}

	double _t;
	SpeciesVector _sv;

	// The fraction of ortho H2 might be important too. Default it here to .25, which
	// corresponds to a ratio of 3 to 1
	double _orthoH2{.75};
} GasStruct;

#endif // CORE_GASSTRUCT_HPP
