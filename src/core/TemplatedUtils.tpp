#include "Error.h"

#include <numeric>

namespace TemplatedUtils
{

/** Return the index permutation required to sort the given data */
template <typename T> std::vector<size_t> argsort(const T* begin, size_t size)
{
	// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes

	// We will sort a vector of indices based on the given range of values
	std::vector<size_t> iv(size);
	std::iota(std::begin(iv), std::end(iv), 0);

	// Sort in-place using a lambda which simply compares the values of the given container
	std::sort(std::begin(iv), std::end(iv),
	          [&](size_t i1, size_t i2) { return begin[i1] < begin[i2]; });

	return iv;
}

/** Return true if value is found in container. */
template <typename T, typename T1> bool contains(const T& value, const T1& container)
{
	return std::find(std::begin(container), std::end(container), value) !=
	       std::end(container);
}

/** Return true if min <= value <= max. */
template <typename T> bool inRange(T value, T min, T max)
{
	return value >= min && value <= max;
}

template <typename T> bool equalWithinTolerance(T value, T reference, T precision)
{
	double minus = reference * (1 - precision);
	double plus = reference * (1 + precision);
	double mn = std::min(minus, plus);
	double mx = std::max(minus, plus);
	return inRange<T>(value, mn, mx);
}

/** Performs a search by bisection on a continuous interval between xMax and xMin. The function
    given as argument should be a decreasing function which indicates in which direction the
    search should proceed, and hence has a zero at the final value of T. The bisection will
    start from xInit, and continues until the final value of x is constrained to an interval of
    width xTolerance. */
template <typename T>
T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax,
                       T xMin)
{
	if (xInit > xMax || xInit < xMin)
		throw "xInit is out of given range for binary search";

	double upperbound = xMax;
	double lowerbound = xMin;
	double current = xInit;

	while (upperbound - lowerbound > xTolerance)
	{
		// Indicates whether the (unknown) final value lies above or below the current
		// one
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

/** Only works with containers which support the insert function */
template <typename T, typename T1> void inline sortedInsert(T elem, T1& container)
{
	container.insert(std::upper_bound(std::begin(container), std::end(container), elem),
	                 elem);
}

/** Evaluate the linear interpolation f(x) of a function f represented by fContainer =
    f(xContainer) */
template <typename T, typename T1, typename T2>
T evaluateLinInterpf(T x, const T1& xContainer, const T2& fContainer)
{
	size_t iRight = index(x, xContainer);
	if (iRight == 0)
		iRight = 1;
	else if (iRight == xContainer.size())
		iRight--;
	size_t iLeft = iRight - 1;
	double val = interpolateLinear(x, xContainer[iLeft], xContainer[iRight],
	                               fContainer[iLeft], fContainer[iRight]);
	return val;
}

template <typename T> T interpolateLinear(T x, T xLeft, T xRight, T fLeft, T fRight)
{
	T weightRight = (x - xLeft) / (xRight - xLeft);
	T val = fLeft + weightRight * (fRight - fLeft);
	return val;
}

template <typename T>
T interpolateRectangular(T x, T y, T xLeft, T xRight, T yLow, T yUp, T fLowerLeft,
                         T fLowerRight, T fUpperLeft, T fUpperRight)
{
	T weightRight = (x - xLeft) / (xRight - xLeft);
	T fLowerI = fLowerLeft + weightRight * (fLowerRight - fLowerLeft);
	T fUpperI = fUpperLeft + weightRight * (fUpperRight - fUpperLeft);

	T weightUpper = (y - yLow) / (yUp - yLow);
	return fLowerI + weightUpper * (fUpperI - fLowerI);
}

/** Trapezoidal integration over two containers. xContainer and yContainer should be the same
    size, and support iteration using std::begin and std::end. Should work with vector and
    valarray for T1 and T2. */
template <typename T, typename T1, typename T2>
T integrate(const T1& xContainer, const T2& yContainer, size_t iMin, size_t iMax)
{
	T answer = 0.0;
	if (xContainer.size() > 1)
	{
		// Start at the right side of the first interval (iMin, iMin + 1)
		auto ix = std::begin(xContainer) + iMin + 1;
		auto iy = std::begin(yContainer) + iMin + 1;
		// One past the last point
		auto ixEnd = std::begin(xContainer) + iMax + 1;
		for (; ix != ixEnd; ix++, iy++)
			answer += 0.5 * (*ix - *(ix - 1)) * (*iy + *(iy - 1));
	}
	else
		answer = yContainer[0];
	return answer;
}

template <typename T, typename T1, typename T2>
T integrate(const T1& xContainer, const T2& yContainer)
{
	return integrate<T, T1, T2>(xContainer, yContainer, size_t{0},
	                            size_t{xContainer.size() - 1});
}

/** Integrates a function over an interval, using an equally spaced number of intermediate
    points. 0 means that only xMin and xMax are used. **/
template <typename T>
T integrateFunction(std::function<T(T)> f, T xMin, T xMax, size_t numPoints)
{
	// 2 extra points |--.--.--| --> three segments
	T xDelta = (xMax - xMin) / static_cast<double>(numPoints + 1);

	std::vector<T> xv(2 + numPoints);
	xv[0] = xMin;
	for (size_t i = 1; i <= numPoints; i++)
		xv[i] = xMin + i * xDelta;
	xv[xv.size() - 1] = xMax;

	std::vector<T> integrandv(xv.size());
	for (size_t i = 0; i < xv.size(); i++)
		integrandv[i] = f(xv[i]);

	return integrate<T, std::vector<T>, std::vector<T>>(xv, integrandv);
}

/** Finds the index of the first element >= val, in the iterable container `container'. If val <
    all the elements, 0 is returned. If val > all the elements, then the index to the last
    element is returned instead. */
template <typename T, typename T1> inline size_t index(T val, const T1& container)
{
	auto iterator = std::lower_bound(std::begin(container), std::end(container), val);
	size_t index = std::distance(std::begin(container), iterator);
	if (index == container.size())
		index--;
	return index;
}

/** Evaluates the polynomial \Sum_{i = 0}^{coeffv.size() - 1} x^i coeffv[i]. */
template <typename T> T evaluatePolynomial(T x, const std::vector<T>& coeffv)
{
	T sum = 0;
	T xToThei = 1;
	for (size_t i = 0; i < coeffv.size(); i++)
	{
		sum += xToThei * coeffv[i];
		xToThei *= x;
	}
	return sum;
}

/** Resamples a function f(x), given on the know points x, to the new points u. This template
    function will work with any combination of types of containers that supports the size() and
    the square bracket dereference operator []. Note that this function has not been tested for
    other types than containers of doubles. The first and second argument provide the known f(x)
    and x values, and should be of equal size. The type of the result can be chosen. The last
    two arguments determine different extrapolation behaviours. -1 means setting everything out
    of the original range to zero. 0 means staying constant once the last point has been passed.
    Any other number means a power law lastvalue * (newpoint / lastpoint), with the argument
    used as an exponent. Stolen from DIRTY/NumUtils.h */
template <typename T, typename T1, typename T2, typename T3>
T linearResample(const T1& function, const T2& knownPoints, const T3& newPoints, int LoEx,
                 int HiEx)
{
	if (function.size() != knownPoints.size())
		Error::runtime("Abscissa and ordinate lengths do not match in interpol");

	T r(newPoints.size());
	double slp, extrslp1 = 0;
	double intcpt, extrint1 = 0;
	size_t j;
	size_t n = knownPoints.size() - 1;
	if (HiEx == -99)
	{
		extrslp1 = (function[n] - function[n - 1]) /
		           (knownPoints[n] - knownPoints[n - 1]);
		extrint1 = function[n] - extrslp1 * knownPoints[n];
	}

	for (size_t i = 0; i < newPoints.size(); ++i)
	{
		if (newPoints[i] < knownPoints[0])
		{ // extrapolate to left.
			switch (LoEx)
			{
			case -1:
				r[i] = 0.;
				break;
			case 0:
				r[i] = function[0];
				break;
			default:
				r[i] = function[0] * pow(newPoints[i] / knownPoints[0], LoEx);
			}
		}
		else if (newPoints[i] > knownPoints[knownPoints.size() - 1])
		{ // extrapolate right.
			switch (HiEx)
			{
			case -1:
				r[i] = 0.;
				break;
			case 0:
				r[i] = function[function.size() - 1];
				break;
			case -99:
				r[i] = extrint1 + extrslp1 * newPoints[i];
				break;
			default:
				r[i] = function[function.size() - 1] *
				       pow(newPoints[i] / knownPoints[knownPoints.size() - 1],
				           HiEx);
				break;
			}
		}
		else
		{
			j = 0;
			while (newPoints[i] > knownPoints[j])
			{
				j++;
			}
			// take care of the case where i=j=0
			// seemed to work w/o this line on 32bit, but not 64bit (even
			// unoptimized)
			if (j == 0)
				j++;
			slp = (function[j] - function[j - 1]) /
			      (knownPoints[j] - knownPoints[j - 1]);
			intcpt = function[j] - slp * knownPoints[j];
			r[i] = slp * newPoints[i] + intcpt;
		}
	}
	return r;
}

} /* namespace TemplatedUtils */
