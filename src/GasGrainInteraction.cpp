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
			double n_d{gsd._grainDensityv[i]};

			// Cross section of the grain
			double sigma_d{gsd._grainSizev[i]};
			sigma_d *= sigma_d / 4.;

			// Formation efficiency epsilon
			double Es;
			double Td; // TODO: graininfo should contain some kind of dust temperature
			double EHp;
			double EHc;
			double a;
			double F;
			double nu_Hc;
			double sqrtEHp_Es = sqrt(EHp - Es);
			double sqrtEHc_Es = sqrt(EHc - Es);

			// 1 / B
			double oneOverB{4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es +
			                8 * sqrt(Constant::PI * Td) / (EHc - Es) *
			                                exp(-2 * a *
			                                    sqrt(2 * Constant::HMASS_CGS /
			                                         Constant::HBAR / Constant::HBAR)) *
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
			total += gsd._grainDensityv[i] * sigma_d * epsilon / oneOverStick;
		}
	}
	return 0.5 * thermalVelocityH;
}
