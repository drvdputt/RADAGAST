#include "GasGrainInteraction.h"
#include "Constants.h"
#include "GrainInfo.h"

double GasGrain::surfaceH2FormationRateCoeff(const GasModule::GrainInfo& g, double thermalVelocityH,
                                             double Tgas)
{
	double total{0};

	for (const auto& gsd : {g._carSizeDist, g._silSizeDist})
	{
		int numSizes = gsd._grainDensityv.size();
		for (int i = 0; i < numSizes; i++)
		{
			// Number density
			double nd{gsd._grainDensityv[i]};

			// Cross section of the grain
			double sigmad{gsd._grainSizev[i]};
			sigmad *= sigmad / 4.;

			// Formation efficiency epsilon
			double Es{gsd._par._es};
			double Td{gsd._tv[i]};
			double EHp{gsd._par._eHp};
			double EHc{gsd._par._eHc};
			double aSqrt{gsd._par._aSqrt};
			double F{gsd._par._f};
			double nu_Hc{gsd._par._nuHc};
			double sqrtEHp_Es = sqrt(EHp - Es);
			double sqrtEHc_Es = sqrt(EHc - Es);

			// 1 / B
			double oneOverB{4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es +
			                8 * sqrt(Constant::PI * Td) / (EHc - Es) *
			                                exp(-2 * aSqrt) *
			                                exp(EHp / Td) * sqrtEHc_Es};
			double onePlusSqrtFrac{1 * sqrtEHc_Es / sqrtEHp_Es};
			double oneOverKsi{1 +
			                  nu_Hc / 2. / F * exp(-1.5 * EHc / Td) * onePlusSqrtFrac *
			                                  onePlusSqrtFrac};
			// eps = (1 + B)^-1 * ksi
			double epsilon{1 / (1 + 1 / oneOverB) / oneOverKsi};

			// Sticking coefficient
			double oneOverStick{1 + 0.04 * sqrt((Tgas + Td) / 100.) +
			                    0.2 * Tgas / 100. + 0.08 * (Tgas * Tgas / 10000.)};

			// factor n_d * sigma_d * epsilon_H2 * S_h
			total += nd * sigmad * epsilon / oneOverStick;
		}
	}
	return 0.5 * thermalVelocityH;
}
