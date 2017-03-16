/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TABLEIMPL_HPP
#define TABLEIMPL_HPP

#include <cstddef>
#include <type_traits>

using std::size_t;

////////////////////////////////////////////////////////////////////

// private compile-time utilities for Table and ArrayTable
namespace Table_Impl
{
    // recursive implementation of "all_off" algorithm
    constexpr bool all() { return true; }
    template <typename... Ts>
    constexpr bool all(bool b, Ts... bs) { return b && all(bs...); }

    // returns true if the parameter pack has the specified number of integer arguments
    template <size_t N, typename... Ts>
    inline constexpr bool isValidArgList()
    {
        return sizeof...(Ts) == N && all(std::is_integral<Ts>::value...);
    }
}

////////////////////////////////////////////////////////////////////

#endif
