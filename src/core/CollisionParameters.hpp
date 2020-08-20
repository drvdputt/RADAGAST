#ifndef CORE_GASSTRUCT_HPP
#define CORE_GASSTRUCT_HPP

#include "EigenAliases.hpp"
#include "SpeciesIndex.hpp"

namespace RADAGAST
{
    /** A set of parameters commonly needed to calculate collisions coefficients. Not every
        function that calculates collision coefficients will need all of this information, but
        using this struct as an argument makes it possible to make those functions virtual. */
    typedef struct CollisionParameters
    {
        CollisionParameters(double t, const SpeciesVector& sv) : _t{t}, _sv{sv} {}
        CollisionParameters(double t, const SpeciesVector& sv, double orthoH2) : _t{t}, _sv{sv}, _orthoH2{orthoH2} {}

        double _t;
        SpeciesVector _sv;

        // The fraction of ortho H2 might be important too. Default it here to .75, which
        // corresponds to the typical ratio of 3 to 1
        double _orthoH2{.75};
    } CollisionParameters;
}
#endif  // CORE_GASSTRUCT_HPP
