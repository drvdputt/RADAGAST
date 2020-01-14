#ifndef CORE_LOOKUPTABLE_HPP
#define CORE_LOOKUPTABLE_HPP

#include "Table.hpp"
#include <string>
#include <vector>

class LookupTable
{
public:
    LookupTable() = default;

    /** Create a lookup table for a set of functions y_i(x). First argument: x-values. Second
        argument: Table containing y_i(xv) on each row (dimension should be (numFunctions, numX). */
    LookupTable(const Array& xv, const Table<2>& yv) : _xv{xv}, _yv{yv} {}

    /** From data file in repo. Needs to be simple column file (see SimpleColumnFile) with the data
        in the first two columns. The resulting lookup table will have one y_i per column (depends
        on numCols given).*/
    LookupTable(const std::string& fname, int numCols, int guessSize);

    /** Do a binary search and evaluate y_i by linearly interpolating for x. */
    double evaluate(int i, double x) const;

private:
    Array _xv;
    Table<2> _yv;
};

#endif  // CORE_LOOKUPTABLE_HPP
