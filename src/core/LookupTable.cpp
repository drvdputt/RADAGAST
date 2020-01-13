#include "LookupTable.hpp"
#include "SimpleColumnFile.hpp"
#include "TemplatedUtils.hpp"

LookupTable::LookupTable(const std::string& fname, int numCols, int guessSize)
{
    if (numCols < 2) Error::runtime("numCols should be >= 2");

    SimpleColumnFile scf(fname);
    scf.read(numCols, guessSize);

    const auto& c0 = scf.column(0);
    _xv = Array(c0.data(), c0.size());

    _yv.resize(numCols - 1, c0.size());
    for (int i = 0; i < _yv.size(0); i++)
    {
        // column zero was already used for x, so start from 1 here
        const auto& ci = scf.column(i + 1);
        for (int j = 0; j < ci.size(); j++) _yv(i, j) = ci[j];
    }
}

double LookupTable::evaluate(int i, double x) const
{
    auto iLeftRight = TemplatedUtils::saneIndexPair(x, _xv);
    return TemplatedUtils::interpolateLinear(x, _xv[iLeftRight[0]], _xv[iLeftRight[1]], _yv(i, iLeftRight[0]),
                                             _yv(i, iLeftRight[1]));
}
