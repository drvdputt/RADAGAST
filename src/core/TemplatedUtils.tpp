#include "Error.hpp"
#include <array>
#include <numeric>
#include <utility>

namespace GasModule
{
    namespace TemplatedUtils
    {
        template<typename T, typename T1> void inline sortedInsert(T elem, T1& container)
        {
            container.insert(std::upper_bound(std::begin(container), std::end(container), elem), elem);
        }

        template<typename T, typename T1> bool contains(const T& value, const T1& container)
        {
            return std::find(std::begin(container), std::end(container), value) != std::end(container);
        }

        template<typename T> bool inRange(T value, T min, T max) { return value >= min && value <= max; }

        template<typename T> bool equalWithinTolerance(T value, T reference, T precision)
        {
            double minus = reference * (1 - precision);
            double plus = reference * (1 + precision);
            double mn = std::min(minus, plus);
            double mx = std::max(minus, plus);
            return inRange<T>(value, mn, mx);
        }

        template<typename T>
        T binaryIntervalSearch(std::function<int(T)> searchDirection, T xInit, T xTolerance, T xMax, T xMin)
        {
            if (xInit > xMax || xInit < xMin) throw "xInit is out of given range for binary search";

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

        template<typename T, typename T1> std::pair<size_t, size_t> saneIndexPair(T x, const T1& xContainer)
        {
            size_t iRight = index(x, xContainer);
            if (iRight == 0)
                iRight = 1;
            else if (iRight == xContainer.size())
                iRight--;
            size_t iLeft = iRight - 1;
            return std::pair<size_t, size_t>(iLeft, iRight);
        }

        template<typename T, typename T1, typename T2>
        T evaluateLinInterpf(T x, const T1& xContainer, const T2& fContainer)
        {
            auto iLeftRight = saneIndexPair(x, xContainer);
            double val = interpolateLinear(x, xContainer[iLeftRight.first], xContainer[iLeftRight.second],
                                           fContainer[iLeftRight.first], fContainer[iLeftRight.second]);
            return val;
        }

        template<typename T> T interpolateLinear(T x, T xLeft, T xRight, T fLeft, T fRight)
        {
            T weightRight = (x - xLeft) / (xRight - xLeft);
            T val = fLeft + weightRight * (fRight - fLeft);
            return val;
        }

        template<typename T>
        T interpolateBilinear(T x, T y, T xLeft, T xRight, T yLow, T yUp, T fLowerLeft, T fLowerRight, T fUpperLeft,
                              T fUpperRight)
        {
            T weightRight = (x - xLeft) / (xRight - xLeft);
            T fLowerI = fLowerLeft + weightRight * (fLowerRight - fLowerLeft);
            T fUpperI = fUpperLeft + weightRight * (fUpperRight - fUpperLeft);

            T weightUpper = (y - yLow) / (yUp - yLow);
            return fLowerI + weightUpper * (fUpperI - fLowerI);
        }

        template<typename T, typename T1, typename T2>
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
                for (; ix != ixEnd; ix++, iy++) answer += 0.5 * (*ix - *(ix - 1)) * (*iy + *(iy - 1));
            }
            else
                answer = yContainer[0];
            return answer;
        }

        template<typename T, typename T1, typename T2> T integrate(const T1& xContainer, const T2& yContainer)
        {
            return integrate<T, T1, T2>(xContainer, yContainer, size_t{0}, size_t{xContainer.size() - 1});
        }

        template<typename T> T integrateFunction(std::function<T(T)> f, T xMin, T xMax, size_t numPoints)
        {
            // 2 extra points |--.--.--| --> three segments
            T xDelta = (xMax - xMin) / static_cast<double>(numPoints + 1);

            std::vector<T> xv(2 + numPoints);
            xv[0] = xMin;
            for (size_t i = 1; i <= numPoints; i++) xv[i] = xMin + i * xDelta;
            xv[xv.size() - 1] = xMax;

            std::vector<T> integrandv(xv.size());
            for (size_t i = 0; i < xv.size(); i++) integrandv[i] = f(xv[i]);

            return integrate<T, std::vector<T>, std::vector<T>>(xv, integrandv);
        }

        template<typename T, typename T1> inline size_t index(T val, const T1& container)
        {
            auto iterator = std::lower_bound(std::begin(container), std::end(container), val);
            size_t index = std::distance(std::begin(container), iterator);
            if (index == container.size()) index--;
            return index;
        }

        template<typename T> T evaluatePolynomial(T x, const std::vector<T>& coeffv)
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
    } /* namespace TemplatedUtils */
}
