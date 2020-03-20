#ifndef CORE_TEMPLATEDUTILS_HPP
#define CORE_TEMPLATEDUTILS_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace GasModule
{
    namespace TemplatedUtils
    {
        /** Insert element into container, keeping it sorted. Only works with containers which
            support the insert function. Should not be used in performance critical code. */
        template<typename T, typename T1> void inline sortedInsert(T elem, T1& container);

        /** Return true if value is found in container. */
        template<typename T, typename T1> bool contains(const T& value, const T1& container);

        /** Return true if min <= value <= max. */
        template<typename T> bool inRange(T value, T min, T max);

        /** Return true if reference * (1 - precision) <= value <= reference * (1 + precision). */
        template<typename T> bool equalWithinTolerance(T value, T reference, T precision);

        /** Performs a search by bisection on a continuous interval between xMax and xMin. The
            function given as argument should be a decreasing function which indicates in which
            direction the search should proceed, and hence has a zero at the final value of T. The
            bisection will start from xInit, and continues until the final value of x is
            constrained to an interval of width xTolerance. */
        template<typename T>
        T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax, T xMin);

        /** Convenience function. Return the index left and right of x, or 0 and 1 if x <
            xContainer[0],q or size-2 and size-1 if x > xContainer[size-1]. */
        template<typename T, typename T1> std::pair<size_t, size_t> saneIndexPair(T x, const T1& xContainer);

        /** Evaluate the linear interpolation f(x) of a function f represented by fContainer =
            f(xContainer) */
        template<typename T, typename T1, typename T2>
        T evaluateLinInterpf(T x, const T1& xContainer, const T2& fContainer);

        template<typename T> T interpolateLinear(T x, T xLeft, T xRight, T fLeft, T fRight);

        template<typename T>
        T interpolateBilinear(T x, T y, T xLeft, T xRight, T yLow, T yUp, T fLowerLeft, T fLowerRight, T fUpperLeft,
                              T fUpperRight);

        /** Trapezoidal integration over two containers. xContainer and yContainer should be the
            same size, and support iteration using std::begin and std::end. Should work with vector
            and valarray for T1 and T2. The last two arguments limit the integration range. The
            first trapezoid is [x[iMin], x[iMin+1]]. The last one is [x[iMax - 1, iMax]]. */
        template<typename T, typename T1, typename T2>
        T integrate(const T1& xContainer, const T2& yContainer, size_t iMin, size_t iMax);

        template<typename T, typename T1, typename T2> T integrate(const T1& xContainer, const T2& yContainer);

        /** Integrates a function over an interval, using an equally spaced number of intermediate
            points. 0 means that only xMin and xMax are used. */
        template<typename T> T integrateFunction(std::function<T(T)> f, T xMin, T xMax);

        /** Finds the index of the first element >= val, in the iterable container `container'. If
            val < all the elements, 0 is returned. If val > all the elements, then the index to the
            last element is returned instead. */
        template<typename T, typename T1> inline size_t index(T val, const T1& container);

        /** Evaluates the polynomial \Sum_{i = 0}^{coeffv.size() - 1} x^i coeffv[i]. */
        template<typename T> T evaluatePolynomial(T x, const std::vector<T>& coeffv);
    } /* namespace TemplatedUtils */
}

#include "TemplatedUtils.tpp"

#endif  // CORE_TEMPLATEDUTILS_HPP
