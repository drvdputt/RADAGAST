#ifndef _TEMPLATEDUTILITIES_H_
#define _TEMPLATEDUTILITIES_H_

#include <algorithm>
#include <functional>
#include <iterator>
#include <valarray>
#include <vector>

namespace TemplatedUtils
{

/* Performs a search by bisection on a continuous interval between xMax and xMin. The function given as
argument should be a decreasing function which indicates in which direction the search should proceed, and
hence has a zero at the final value of T. The bisection will start from xInit, and continues until the final
value of x is constrained to an interval of width xTolerance. */
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

/* Only works with containers which support the insert function */
template<typename T, typename T1>
void inline sortedInsert(T elem, T1& container)
{
	container.insert(std::upper_bound(std::begin(container), std::end(container), elem), elem);
}

template<typename T>
T interpolateRectangular(T x, T y, T xLeft, T xRight, T yLow, T yUp, T fLowerLeft, T fLowerRight,
		T fUpperLeft, T fUpperRight)
{
	T weightRight = (x - xLeft) / (xLeft - xRight);
	T fLowerI = fLowerLeft + weightRight * (fLowerRight - fLowerLeft);
	T fUpperI = fUpperLeft + weightRight * (fUpperRight - fUpperLeft);

	T weightUpper = (y - yLow) / (yUp - yLow);
	return fLowerI + weightUpper * (fUpperI - fLowerI);
}

/* Trapezoidal integration over two containers. xContainer and yContainer should be the same size, and support
iteration using std::begin and std::end. Should work with vector and valarray for T1 and T2. */
template<typename T, typename T1, typename T2>
T integrate(const T1& xContainer, const T2& yContainer)
{
	T answer = 0.0;
	if (xContainer.size() > 1)
	{
		auto iy = std::begin(yContainer) + 1;
		for (auto ix = std::begin(xContainer) + 1; ix != std::end(xContainer); ix++, iy++)
		{
			answer += 0.5 * (*ix - *(ix - 1)) * (*iy + *(iy - 1));
		}
	}
	else
	{
		answer = yContainer[0];
	}
	return answer;
}

/* Finds the index of the first element >= val, in the iterable container `container'. */
template<typename T, typename T1> inline size_t index(T val, const T1& container)
{
	auto idx = std::lower_bound(std::begin(container), std::end(container), val);
	return std::distance(std::begin(container), idx);
}

} /* namespace TemplatedUtils */

#endif /* _TEMPLATEDUTILITIES_H_ */
