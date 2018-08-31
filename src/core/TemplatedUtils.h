#ifndef _TEMPLATEDUTILITIES_H_
#define _TEMPLATEDUTILITIES_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace TemplatedUtils
{
template <typename T> bool inRange(T value, T min, T max);

template <typename T> bool equalWithinTolerance(T value, T reference, T precision);

template <typename T>
T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax,
                       T xMin);

template <typename T, typename T1> void inline sortedInsert(T elem, T1& container);

template <typename T, typename T1, typename T2>
T evaluateLinInterpf(T x, const T1& xContainer, const T2& fContainer);

template <typename T> T interpolateLinear(T x, T xLeft, T xRight, T fLeft, T fRight);

template <typename T>
T interpolateRectangular(T x, T y, T xLeft, T xRight, T yLow, T yUp, T fLowerLeft,
                         T fLowerRight, T fUpperLeft, T fUpperRight);

template <typename T, typename T1, typename T2>
T integrate(const T1& xContainer, const T2& yContainer, size_t iMin, size_t iMax);

template <typename T, typename T1, typename T2>
T integrate(const T1& xContainer, const T2& yContainer);

template <typename T> T integrateFunction(std::function<T(T)> f, T xMin, T xMax);

template <typename T, typename T1> inline size_t index(T val, const T1& container);

template <typename T> T evaluatePolynomial(T x, const std::vector<T>& coeffv);

template <typename T, typename T1, typename T2, typename T3>
T linearResample(const T1& function, const T2& knownPoints, const T3& newPoints, int LoEx,
                 int HiEx);
} /* namespace TemplatedUtils */

#include "TemplatedUtils.tpp"

#endif /* _TEMPLATEDUTILITIES_H_ */
