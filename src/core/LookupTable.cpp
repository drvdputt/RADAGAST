#include "LookupTable.hpp"
#include "SimpleColumnFile.hpp"
#include "TemplatedUtils.hpp"

namespace RADAGAST
{
    LookupTable::LookupTable(const std::string& fname, int numCols, int guessSize)
    {
        if (numCols < 2) Error::runtime("numCols should be >= 2");

        InColumnFile input(fname);
        input.read(numCols, guessSize);

        const auto& c0 = input.column(0);
        _xv = Array(c0.data(), c0.size());

        _yv.resize(numCols - 1, c0.size());
        for (int i = 0; i < _yv.size(0); i++)
        {
            // column zero was already used for x, so start from 1 here
            const auto& ci = input.column(i + 1);
            for (int j = 0; j < ci.size(); j++) _yv(i, j) = ci[j];
        }
    }

    double LookupTable::evaluate(int i, double x, bool extrapolate) const
    {
        if (!extrapolate && (x < _xv[0] || _xv[_xv.size() - 1] < x)) return 0.;
        auto iLeftRight = TemplatedUtils::saneIndexPair(x, _xv);
        return TemplatedUtils::interpolateLinear(x, _xv[iLeftRight.first], _xv[iLeftRight.second],
                                                 _yv(i, iLeftRight.first), _yv(i, iLeftRight.second));
    }
}
