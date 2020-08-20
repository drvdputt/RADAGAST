#ifndef CORE_LOOKUPTABLE_HPP
#define CORE_LOOKUPTABLE_HPP

#include "Table.hpp"
#include <string>
#include <vector>

namespace RADAGAST
{
    class LookupTable
    {
    public:
        LookupTable() = default;

        /** Create a lookup table for a set of functions y_i(x). First argument: x-values. Second
            argument: Table containing y_i(xv) on each row (dimension should be (numFunctions,
            numX). */
        LookupTable(const Array& xv, const Table<2>& yv) : _xv{xv}, _yv{yv} {}

        /** Create a lookup table from a data file in this repo. Needs to be simple column file
            (see SimpleColumnFile) with at least two columns (and numCols should be >= 2).
            Column 0 will be used as for xv, and the other (numCols - 1) column will get an y_i
            entry (with i = 0 for column 1 of the file).*/
        LookupTable(const std::string& fname, int numCols, int guessSize);

        /** Do a binary search and evaluate y_i by linearly interpolating for x. When the given
            x is out of the bound of the data, return zero if extrapolate is false. Extrapolate
            linearly using the last two points if extrapolate is true. */
        double evaluate(int i, double x, bool extrapolate = false) const;

    private:
        Array _xv;
        Table<2> _yv;
    };
}
#endif  // CORE_LOOKUPTABLE_HPP
