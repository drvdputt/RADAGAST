#include "GasGrainInteraction.h"
#include "Constants.h"
#include "Error.h"
#include "GrainInterface.h"
#include "GrainType.h"

double GasGrain::surfaceH2FormationRateCoeff(const GasModule::GrainInterface& gInterface,
                                             double Tgas)

{
	// See Rollig et al. (2013) appendix C + erratum of 2002 Cazaux and Tielens paper
	double total{0};

	double Tgas_100 = Tgas / 100.;

	size_t numPop = gInterface.numPopulations();
	for (size_t i = 0; i < numPop; i++)
	{
		const GasModule::GrainInterface::Population* pop = gInterface.population(i);
		const GrainType* grainType = pop->type();
		const GasModule::SfcInteractionPar& surfaceParams =
		                grainType->sfcInteractionPar();

		double Es{surfaceParams._es};
		double EHp{surfaceParams._eHp};
		double EHc{surfaceParams._eHc};
		double aSqrt{surfaceParams._aSqrt};
		double F{surfaceParams._f};
		double nu_Hc{surfaceParams._nuHc};

		double EHc_Es = EHc - Es;
		double sqrtEHp_Es = sqrt(EHp - Es);
		double sqrtEHc_Es = sqrt(EHc_Es);
		double sqrtEHc_Ehp = sqrt(EHc - EHp);
		double onePlusSqrtFrac = 1. + sqrtEHc_Es / sqrtEHp_Es;

		size_t numSizes = pop->numSizes();
		for (size_t iSize = 0; iSize < numSizes; iSize++)
		{
			// Number density
			double nd{pop->densityv()[iSize]};
			double Td{pop->temperaturev()[iSize]};

			// Cross section of the grain. This needs to be an integrated quantity
			// over the bin, but lets approximate for now. TODO: integrated cross
			// section for grains.
			double sigmad{pop->sizev()[iSize]};
			sigmad *= sigmad * Constant::PI;

			// 1 / B
			double beta_alpha = 1. / (4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es +
			                          8. * sqrt(Constant::PI * Td) *
			                                          exp(-2. * aSqrt + EHp / Td) *
			                                          sqrtEHc_Ehp / EHc_Es);

			double xi = 1. / (1. + nu_Hc * exp(-1.5 * EHc / Td) * onePlusSqrtFrac *
			                                       onePlusSqrtFrac / 2. / F);
			// The minus sign in the exponential is not there in Rollig 2013! The
			// erratum for Cazaux and Tielens (2002) has the correct version of
			// formula 16 of CT02).

			// eps = (1 + B)^-1 * ksi
			double epsilon = xi / (1. + beta_alpha);

			// Sticking coefficient (same paper, equation 20; originally from
			// Hollenbach and Mckee (!979), equation 3.7, or Burke and Hollenbach
			// (1979)). Annoyingly, Rollig et al states 0.04 instead of 0.4 for the
			// first coefficient.
			double S = 1. / (1. + 0.4 * sqrt(Tgas_100 + Td / 100.) +
			                 0.2 * Tgas_100 + 0.08 * Tgas_100 * Tgas_100);

			// factor n_d * sigma_d * epsilon_H2 * S_h
			double contribution = nd * sigmad * epsilon * S;
			if (std::isfinite(contribution))
				total += contribution;
		}
	}
	// Use mean particle speed of Maxwell distribution, not RMS as suggested by comment in
	// Cloudy source code.
	double thermalVelocityH = sqrt(8. * Constant::BOLTZMAN / Constant::PI * Tgas /
	                               Constant::HMASS_CGS);
	return 0.5 * thermalVelocityH * total;
}

double GasGrain::simpleGasGrainCool(double Tdust, double Tgas, double nH, double nH2)
{
	double lambda_gd = 0;

	double T_factor = sqrt(Tgas) * (Tgas - Tdust);

	// H contribution
	double Cgd = 3.8e-33; // erg s-1 cm3 K
	double n2 = nH * nH;
	lambda_gd += Cgd * T_factor * n2;

	// H2 contribution
	Cgd = 1.0e-33; // erg s-1 cm3 K
	n2 = nH2 * nH2;
	lambda_gd += Cgd * T_factor * n2;

	return lambda_gd;
}
