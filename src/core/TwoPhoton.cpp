#include "TwoPhoton.hpp"
#include "Constants.hpp"

namespace RADAGAST
{
    Array TwoPhoton::emissivityv(const Array& eFrequencyv, double n2s)
    {
        // 1984-Nussbaumer, constant factor in eq 3 and energy difference between 2s and 1s.
        // While this data is also contained in the HData, pasting it here is easier. The number
        // (in cm-1) comes from h_1.elvlc (file from CHIANTI).
        const double nu0 = 82258.956 * Constant::LIGHT;
        const double constFactor = Constant::PLANCK / Constant::FPI * n2s;

        // Parameters for eq 2
        const double C = 202.0;  // s-1
        const double alpha = .88;
        const double beta = 1.53;
        const double gam = .8;

        Array result(eFrequencyv.size());
        for (size_t iFreq = 0; eFrequencyv[iFreq] < nu0 && iFreq < eFrequencyv.size(); iFreq++)
        {
            double y = eFrequencyv[iFreq] / nu0;
            double y1miny = y * (1 - y);
            double pow4y1miny_gam = pow(4 * y1miny, gam);
            double Py = C * (y1miny * (1 - pow4y1miny_gam) + alpha * pow(y1miny, beta) * pow4y1miny_gam);
            result[iFreq] = constFactor * y * Py;
        }
        return result;
    }
}
