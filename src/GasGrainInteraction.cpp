#include "GasGrainInteraction.h"
#include "Constants.h"

#include "GrainInterface.h"
#include "GrainTypeProperties.h"

double GasGrain::surfaceH2FormationRateCoeff(const GasModule::GrainInterface& gInterface,
                                             double Tgas)
{
	double total{0};

	double thermalVelocityH{sqrt(Constant::BOLTZMAN * Tgas / Constant::HMASS_CGS)};

	// Using auto here to avoid having to specify all the namespaces
	const auto& grainPopv = gInterface.populationv();
	for (const auto& grainPop : grainPopv)
	{
		const auto& surfaceParams = GrainTypeProperties::sfcInteractionPar(grainPop._type);
		int numSizes = grainPop._sizev.size();
		for (int iSize = 0; iSize < numSizes; iSize++)
		{
			// Number density
			double nd{grainPop._densityv[iSize]};
			double Td{grainPop._temperaturev[iSize]};

			// Cross section of the grain
			double sigmad{grainPop._sizev[iSize]};
			sigmad *= sigmad / 4.;

			// Formation efficiency epsilon
			double Es{surfaceParams._es};
			double EHp{surfaceParams._eHp};
			double EHc{surfaceParams._eHc};
			double aSqrt{surfaceParams._aSqrt};
			double F{surfaceParams._f};
			double nu_Hc{surfaceParams._nuHc};
			double sqrtEHp_Es = sqrt(EHp - Es);
			double sqrtEHc_Es = sqrt(EHc - Es);

			// 1 / B
			double oneOverB{4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es +
			                8 * sqrt(Constant::PI * Td) / (EHc - Es) * exp(-2 * aSqrt) *
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
	return 0.5 * thermalVelocityH * total;
}
