#include "GrainH2Formation.hpp"
#include "Constants.hpp"
#include "Functions.hpp"

namespace GasModule
{
    double GrainH2Formation::surfaceH2FormationRateCoeff(const Array& sizev, const Array& temperaturev,
                                                         const Array& densityv, double Tgas) const
    {
        // See Rollig et al. (2013) appendix C + erratum of 2002 Cazaux and Tielens paper
        const Array& coeffPerGrainPerHPerSizev = surfaceH2FormationRateCoeffPerSize(sizev, temperaturev, Tgas);
        double total = (densityv * coeffPerGrainPerHPerSizev).sum();
        return total;
    }

    Array GrainH2Formation::surfaceH2FormationRateCoeffPerSize(const Array& sizev, const Array& temperaturev,
                                                               double Tgas) const
    {
        size_t numSizes = sizev.size();
        Array formationPerGrainPerHPerSizev(numSizes);

        double Es{_sfcInteractionPar._es};
        double EHp{_sfcInteractionPar._eHp};
        double EHc{_sfcInteractionPar._eHc};
        double aSqrt{_sfcInteractionPar._aSqrt};
        double F{_sfcInteractionPar._f};
        double nu_Hc{_sfcInteractionPar._nuHc};

        double EHc_Es = EHc - Es;
        double sqrtEHp_Es = sqrt(EHp - Es);
        double sqrtEHc_Es = sqrt(EHc_Es);
        double sqrtEHc_Ehp = sqrt(EHc - EHp);
        double onePlusSqrtFrac = 1. + sqrtEHc_Es / sqrtEHp_Es;

        double Tgas_100 = Tgas / 100.;
        double vH = Functions::meanThermalVelocity(Tgas, Constant::HMASS);

        for (size_t i = 0; i < numSizes; i++)
        {
            // Cross section of the grain. This actually needs to be average(a^2) over the grain
            // bin, and not average(a)^2, but lets approximate with the latter for now.
            double Td{temperaturev[i]};
            double sigmad{sizev[i]};
            sigmad *= sigmad * Constant::PI;

            // 1 / B
            double beta_alpha = 1.
                                / (4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es
                                   + 8. * sqrt(Constant::PI * Td) * exp(-2. * aSqrt + EHp / Td) * sqrtEHc_Ehp / EHc_Es);

            double xi = 1. / (1. + nu_Hc * exp(-1.5 * EHc / Td) * onePlusSqrtFrac * onePlusSqrtFrac / 2. / F);

            // The minus sign in the exponential is not there in Rollig 2013! The erratum for
            // Cazaux and Tielens (2002) has the correct version of formula 16 of CT02).

            // eps = (1 + B)^-1 * ksi
            double epsilon = xi / (1. + beta_alpha);

            // Sticking coefficient (same paper, equation 20; originally from Hollenbach and Mckee
            // (1979), equation 3.7, or Burke and Hollenbach (1979)). Annoyingly, Rollig et al
            // states 0.04 instead of 0.4 for the first coefficient.
            double S = 1. / (1. + 0.4 * sqrt(Tgas_100 + Td / 100.) + 0.2 * Tgas_100 + 0.08 * Tgas_100 * Tgas_100);

            // sigma_d * epsilon_H2 * S_h. Needs to be multiplied with the grain number density
            // later
            double product = 0.5 * vH * sigmad * epsilon * S;
            if (std::isfinite(product)) formationPerGrainPerHPerSizev[i] = product;
            // Else it will stay 0
        }
        return formationPerGrainPerHPerSizev;
    }

    Array GrainH2Formation::surfaceH2FormationHeatPerSize(const Array& sizev, const Array& temperaturev, double Tgas,
                                                          double nH) const
    {
        const Array& coeffPerGrainPerHPerSizev = surfaceH2FormationRateCoeffPerSize(sizev, temperaturev, Tgas);
        return coeffPerGrainPerHPerSizev * nH * _heatPerH2;
    }
}
