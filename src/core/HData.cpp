#include "HData.hpp"
#include "Constants.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "RecombinationRate.hpp"

namespace RADAGAST
{
    HData::HData() : LevelCoefficients(Constant::HMASS) {}

    EVector HData::recombinationRatev(double T) const
    {
        EVector result{EVector::Zero(numLv())};
        // Now loop over all levels, and add the correct recombination coefficient. If a level
        // is collapsed, the alpha will be added to the same level multiple times.
        for (int n = 1; n <= nMax(); n++)
        {
            for (int l = 0; l < n; l++)
            {
                double maovalue = _rr.alpha(n, l, T);
                result[index(n, l)] += maovalue;
            }
        }
        if (Options::hlevels_topoff)
        {
            // The recombination coefficients should sum to the total one (the one used for the
            // ionization balance calculation). Therefore, we calculate here how much we still
            // need...
            double total_rr = Ionization::recombinationRateCoeff(T);
            double topoff = total_rr - result.sum();
            if (topoff > 0)
            {
                // ...and add it to the top n-level
                for (int l = 0; l < nMax(); l++)
                {
                    // These weights should sum to 1
                    double weight = (2 * l + 1) / nMax() / nMax();
                    result[index(nMax(), l)] += topoff * weight;
                }
            }
        }
        return result;
    }
}
