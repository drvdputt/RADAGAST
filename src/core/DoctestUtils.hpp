#ifndef CORE_DOCTESTUTILS_HPP
#define CORE_DOCTESTUTILS_HPP

#include "EigenAliases.hpp"
#include "TemplatedUtils.hpp"
#include "doctest.h"
#include <sstream>

namespace RADAGAST
{
    /** Here are some wrappers around some of the doctest assertion macros. Frequently used checks
        are defined here, with some standardized, readable failure messages. */
    namespace DoctestUtils
    {
        // Option that makes some of my tests dump data to files, to do some manual
        // checking/plotting of results.
        const bool allowFileOutput = false;

        inline void checkRange(std::string quantityName, double value, double min, double max, bool warn = false)
        {
            bool inRange = TemplatedUtils::inRange(value, min, max);
            if (!inRange)
            {
                std::stringstream ss;
                ss << quantityName << " = " << value << ". Should be between " << min << " and " << max;
                std::string message = ss.str();
                if (warn)
                    WARN_MESSAGE(inRange, message);
                else
                    CHECK_MESSAGE(inRange, message);
            }
        }

        inline void checkTolerance(std::string quantityName, double value, double reference, double precision,
                                   bool warn = false)
        {
            bool withinTolerance = TemplatedUtils::equalWithinTolerance<double>(value, reference, precision);
            if (!withinTolerance)
            {
                std::stringstream ss;
                ss << quantityName << " " << value << " is not within " << precision << " precision of reference value "
                   << reference;
                std::string message = ss.str();
                if (warn)
                    WARN_MESSAGE(withinTolerance, message);
                else
                    CHECK_MESSAGE(withinTolerance, message);
            }
        }

        template<typename Derived>
        void compareMatrices(const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<Derived>& b, double tolerance)
        {
            // Calculate relative difference element-wise
            EMatrix relDiff = (a - b).array() / b.array();
            for (int i = 0; i < relDiff.size(); i++)
            {
                // Take out the nans (due to divide by zero)
                auto* pointer = relDiff.data() + i;
                auto value = *pointer;
                *pointer = std::isfinite(value) ? value : 0;
            }
            bool withinTolerance = (relDiff.cwiseAbs().array() < tolerance).all();
            if (!withinTolerance)
            {
                std::stringstream ss;
                ss << "matrix " << a << " is not within " << tolerance << " precision of reference value " << b;
                std::string message = ss.str();
                CHECK_MESSAGE(withinTolerance, message);
            }
        }
    }  // namespace DoctestUtils
}
#endif  // CORE_DOCTESTUTILS_HPP
