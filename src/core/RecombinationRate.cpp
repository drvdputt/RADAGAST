#include "RecombinationRate.hpp"
#include "Constants.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "TemplatedUtils.hpp"
#include <string>

namespace
{
    const std::map<char, int> spdf_to_l = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3},
                                           {'G', 4}, {'H', 5}, {'I', 6}, {'J', 7}};
    int NUMROWSMAO = 272; // 256 data rows + 16 unused (s states need only one row)
}

namespace GasModule
{
    HydrogenADF48::HydrogenADF48() { readADF48File("dat/adf48-cut"); }

    HydrogenADF48::~HydrogenADF48() = default;

    void HydrogenADF48::readADF48File(const std::string& path)
    {
        auto ifs = IOTools::ifstreamRepoFile(path);
        std::string line;

        // indexing
        getline(ifs, line);
        int index = 0;
        while (getline(ifs, line))
        {
            if (line.empty()) break;
            int i, g;
            std::string nl_string;
            std::istringstream(line) >> i >> g >> nl_string;
            int n = nl_string[0] - '0';
            int l = spdf_to_l.at(nl_string[1]);
            _nlToIndexm.insert({{n, l}, index});
            index++;
        }
        // temperature
        getline(ifs, line);
        getline(ifs, line);
        auto tv = IOTools::allNumbersFromNextLine(line);
        _temperaturev = Array(tv.data(), tv.size());

        // recombination coefficient
        getline(ifs, line);
        while (getline(ifs, line))
        {
            std::vector<double> alphav = IOTools::allNumbersFromNextLine(line);
            Error::equalCheck("Number of columns in adf48 file", alphav.size(), tv.size());
            _alphavv.emplace_back(alphav);
        }
    }

    double HydrogenADF48::alpha(int n, int l, double T) const
    {
        // Grab the right line
        const std::vector<double>& alphav = _alphavv[_nlToIndexm.at({n, l})];

        // Interpolate on the temperature
        return TemplatedUtils::evaluateLinInterpf(T, _temperaturev, alphav);
    }

    Mao2016RecombinationRate::Mao2016RecombinationRate()
    {
        // data goes up to n = 16, so the number of lines (nlj combinations) is 256, and the
        // number of nl combinations is 128 (2 different j per nl).
        _fitCoefficients.resize(NUMROWSMAO, 7);

        auto ifs = IOTools::ifstreamRepoFile("dat/rr_spex3(1).dat");
        std::string line;
        int linenr = 0;

        while (getline(ifs, line))
        {
            linenr += 1;
            if (line.front() == '#') continue;

            // read the values as listed in the readme (I pasted the readme as a comment at the
            // top of the file)

            // isosequence and atomic number
            int s, z;
            // fit parameters and unused numbers
            std::array<double, 7> a0_b0_c0_a1_b1_a2_b2;
            double unused;
            // quantum numbers
            int n, sp, l;
            double j;

            std::istringstream iss(line);
            iss >> s >> z;
            for (int i = 0; i < 7; i++) iss >> a0_b0_c0_a1_b1_a2_b2[i];
            iss >> unused >> unused >> unused >> n >> sp >> l >> j;

            // j is either l - 0.5 or l + 0.5, except for s states. Use the highest index for
            // the second state.
            int tableIndex = 2 * nlToIndex(n, l);
            if (l != 0 && j > l) tableIndex += 1;

            // copy the fit parameters
            for (int i = 0; i < 7; i++) _fitCoefficients(tableIndex, i) = a0_b0_c0_a1_b1_a2_b2[i];
        }
    }

    Mao2016RecombinationRate::~Mao2016RecombinationRate() = default;

    double Mao2016RecombinationRate::alpha(int n, int l, double T) const
    {
        double T_eV = T * Constant::BOLTZMAN / Constant::EV;

        // for both j
        double total = 0;
        int numj = l == 0 ? 1 : 2; // s states can have only one value for j
        for (int offset = 0; offset < numj; offset++)
        {
            // retrieve the parameters
            int i = 2 * nlToIndex(n, l) + offset;
            if (i > NUMROWSMAO - 1) Error::runtime("Out of range for Mao et al. (2016) recombination data");
            double a0 = _fitCoefficients(i, 0);
            double b0 = _fitCoefficients(i, 1);
            double c0 = _fitCoefficients(i, 2);
            double a1 = _fitCoefficients(i, 3);
            double b1 = _fitCoefficients(i, 4);
            double a2 = _fitCoefficients(i, 5);
            double b2 = _fitCoefficients(i, 6);

            // apply the fitting formula (cm3 s-1)
            total += 1e-10 * a0 * std::pow(T_eV, -b0 - c0 * std::log(T_eV)) * (1 + a2 * std::pow(T_eV, -b2))
                     / (1 + a1 * std::pow(T_eV, -b1));
        }
        return total;
    }
}
