#include "LookupTable.hpp"
#include "SimpleColumnFile.hpp"
#include "TemplatedUtils.hpp"

LookupTable::LookupTable(const std::string& fname, int numCols, int guessSize)
{
	SimpleColumnFile scf(fname);
	scf.read(numCols, guessSize);

	const auto& c0 = scf.column(0);
	_xv = Array(c0.data(), c0.size());

	_yv.resize(numCols, c0.size());
	for (int i = 1; i < numCols; i++)
	{
		const auto& ci = scf.column(i);
		for (int j = 0; j < ci.size(); j++)
			_yv(i, j) = ci[j];
	}
}

double LookupTable::evaluate(int i, double x) const
{
	auto iLeftRight = TemplatedUtils::saneIndexPair(x, _xv);
	return TemplatedUtils::interpolateLinear(x, _xv[iLeftRight[0]], _xv[iLeftRight[1]],
	                                         _yv(i, iLeftRight[0]), _yv(i, iLeftRight[1]));
}
