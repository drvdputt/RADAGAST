#ifndef CORE_RECOMBINATIONRATE_HPP
#define CORE_RECOMBINATIONRATE_HPP

#include "Array.hpp"
#include "Table.hpp"
#include <array>
#include <map>
#include <string>
#include <vector>

namespace GasModule
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
}
#endif  // CORE_RECOMBINATIONRATE_HPP
