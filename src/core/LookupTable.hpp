#ifndef CORE_LOOKUPTABLE_HPP
#define CORE_LOOKUPTABLE_HPP

#include <vector>

class LookupTable
{
public:
	/** By giving two lists of values */
	LookupTable(const std::vector<double>& xv, const std::vector<double>& yv)
	                : _xv{xv}, _yv{yv}
	{
	}

	/** From data file in repo. Needs to be simple column file (see SimpleColumnFile) with
	    the data in the first two columns. */
	LookupTable(const std::string& fname, int guessSize);

	/** Do a binary search and evaluate y by linearly interpolating for x. */
	double evaluate(double x) const;

private:
	std::vector<double> _xv;
	std::vector<double> _yv;
};

#endif // CORE_LOOKUPTABLE_HPP
