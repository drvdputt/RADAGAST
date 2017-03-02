#ifndef _TEMPLATEDUTILITIES_H_
#define _TEMPLATEDUTILITIES_H_

namespace TemplatedUtils
{

// Performs a search by bisection on a continuous interval between xMax and xMin. The function given as argument
// should be a decreasing function which indicates in which direction the search should proceed, and hence has
// a zero at the final value of T. The bisection will start from xInit, and continues until the final value of x
// is constrained to an interval of width xTolerance.
template<typename T>
T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax, T xMin)
{
	if (xInit > xMax || xInit < xMin)
		throw "xInit is out of given range for binary search";

	double upperbound = xMax;
	double lowerbound = xMin;
	double current = xInit;

	while (upperbound - lowerbound > xTolerance)
	{
		// Indicates whether the (unknown) final value lies above or below the current one
		double iCompare = searchDirection(current);

		if (iCompare > 0)
			lowerbound = current;
		else if (iCompare < 0)
			upperbound = current;
		else
			return current;

		current = (lowerbound + upperbound) / 2.;
	}
	return current;
}

template<typename T>
void sorted_insert(vector<T>& vec, T elem)
{
	vec.insert(upper_bound(vec.begin(), vec.end(), elem), elem);
}

}

#endif /* _TEMPLATEDUTILITIES_H_ */
