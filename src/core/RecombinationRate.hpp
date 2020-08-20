#ifndef CORE_RECOMBINATIONRATE_HPP
#define CORE_RECOMBINATIONRATE_HPP

#include "Array.hpp"
#include "Table.hpp"
#include <array>
#include <map>
#include <string>
#include <vector>

namespace RADAGAST
{
    /** Abstract class representing the recombination rate of H. */
    class RecombinationRate
    {
    public:
        RecombinationRate() = default;
        virtual ~RecombinationRate() = default;

        /** Calculate the recombination rate to level n,l at temperature T */
        virtual double alpha(int n, int l, double T) const = 0;
    };

    /** Partial radiative recombination rates read from the ADF48 file for Hydrogen. Because
        redistribution of this file is not allowed, the user has to download it first to use
        this class. The data are tabulated for different temperatures, go up to n = 8, and are
        resolved on l. WARNING: currently broken TODO: make sure that this class works with the
        file as downloaded straight from openADAS. The current code only works with a version
        that I edited by hand. */
    class HydrogenADF48 : public RecombinationRate
    {
    public:
        HydrogenADF48();
        ~HydrogenADF48();

    private:
        /** Load the recombination data from a modified ADF48 file (from openADAS). */
        void readADF48File(const std::string& path);

    public:
        /** Interpolate the recombination data of the ADF48 file for a certain temperature. */
        double alpha(int n, int l, double T) const override;

    private:
        std::map<std::array<int, 2>, size_t> _nlToIndexm;
        Array _temperaturev;
        std::vector<std::vector<double>> _alphavv;
    };

    /** Partial radiative recombination rates from a research note by Mao and Kaastra (2016). A
    fitting formula similar to that for the total recombination rate from Verner (1996) is used,
    and the coefficients are provided in a file with data up to n = 16. */
    class Mao2016RecombinationRate : RecombinationRate
    {
    public:
        /** Load the fitting coefficients from the file rr_spex3(1).dat */
        Mao2016RecombinationRate();
        ~Mao2016RecombinationRate();

        /** Apply the fitting formula using the coefficient for the given level */
        double alpha(int n, int l, double T) const override;

    private:
        /** Mapping from nl to an index can be done using a simple formula. The number of levels
            up to and including a certain n is n(n+1) / 2. So an easy formula for a unique index
            is (n - 1)n / 2 + l. Internally, each of these indices i refers to a pair of rows in
            the fit coefficient table, located at 2i and 2i+1, pointing to the j = l - 0.5 and j
            = l + 0.5 levels respectively. The maximum index is (16 + 1) * 16 / 2 - 1 = 135, and
            the fit coefficient table will have 136 * 2 = 272 rowls. This is 16 more than the
            256 data points, because s states only have one j. So for each of the 16 different
            n, there is one data less than the number obtained by doubling the number of nl
            combinations. */
        int nlToIndex(int n, int l) const { return ((n - 1) * n) / 2 + l; }

        // a0 b0 c0 a1 b1 a2 b2 for every n,l (indexed according to the above map). 256 rows in
        // total, as the number of nl combinations is 128, and there are two entries for each
        // because of the 2j+1 quantum number.
        Table<2> _fitCoefficients;
    };
}
#endif  // CORE_RECOMBINATIONRATE_HPP
