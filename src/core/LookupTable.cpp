#include "LookupTable.hpp"
#include "SimpleColumnFile.hpp"
#include "TemplatedUtils.hpp"

LookupTable::LookupTable(const std::string& fname, int guessSize)
{
	SimpleColumnFile scf(fname);
	scf.read(2, guessSize);
	_xv = scf.column(0);
	_yv = scf.column(1);
}

double LookupTable::evaluate(double x) const
{
	return TemplatedUtils::evaluateLinInterpf(x, _xv, _yv);
}
