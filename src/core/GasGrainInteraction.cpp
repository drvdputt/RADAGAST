#include "GasGrainInteraction.hpp"
#include "Error.hpp"
#include "GrainType.hpp"

Array GasGrain::surfaceH2FormationRateCoeffPerSize(
                const GasModule::GrainInterface::Population& pop, double Tgas)
{
	size_t numSizes = pop.numSizes();
	Array formationPerGrainPerHPerSizev(numSizes);

	const GrainType* grainType = pop.type();
	const GasModule::SfcInteractionPar& surfaceParams = grainType->sfcInteractionPar();

	if (!surfaceParams._valid)
	{
		// If this struct was created using the default constructor, return zeros.
		// (Typical use case: a new grain type (implemented later) with unknown
		// properties for H2 formation.)
		return formationPerGrainPerHPerSizev;
	}

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

	double Tgas_100 = Tgas / 100.;

	// Use mean particle speed of Maxwell distribution, not RMS as suggested by comment in
	// Cloudy source code.
	double thermalVelocityH = sqrt(8. * Constant::BOLTZMAN / Constant::PI * Tgas /
	                               Constant::HMASS_CGS);

	for (size_t iSize = 0; iSize < numSizes; iSize++)
	{
		double Td{pop.temperaturev()[iSize]};

		// Cross section of the grain. This needs to be an integrated quantity over the
		// bin, but lets approximate for now. TODO: integrated cross section for grains.
		double sigmad{pop.sizev()[iSize]};
		sigmad *= sigmad * Constant::PI;

		// 1 / B
		double beta_alpha = 1. / (4. * exp(Es / Td) * sqrtEHp_Es / sqrtEHc_Es +
		                          8. * sqrt(Constant::PI * Td) *
		                                          exp(-2. * aSqrt + EHp / Td) *
		                                          sqrtEHc_Ehp / EHc_Es);

		double xi = 1. / (1. + nu_Hc * exp(-1.5 * EHc / Td) * onePlusSqrtFrac *
		                                       onePlusSqrtFrac / 2. / F);

		// The minus sign in the exponential is not there in Rollig 2013! The erratum
		// for Cazaux and Tielens (2002) has the correct version of formula 16 of CT02).

		// eps = (1 + B)^-1 * ksi
		double epsilon = xi / (1. + beta_alpha);

		// Sticking coefficient (same paper, equation 20; originally from Hollenbach and
		// Mckee (1979), equation 3.7, or Burke and Hollenbach (1979)). Annoyingly,
		// Rollig et al states 0.04 instead of 0.4 for the first coefficient.
		double S = 1. / (1. + 0.4 * sqrt(Tgas_100 + Td / 100.) + 0.2 * Tgas_100 +
		                 0.08 * Tgas_100 * Tgas_100);

		// sigma_d * epsilon_H2 * S_h. Needs to be multiplied with the grain number density later
		double product = 0.5 * thermalVelocityH * sigmad * epsilon * S;
		if (std::isfinite(product))
			formationPerGrainPerHPerSizev[iSize] = product;
	}
	return formationPerGrainPerHPerSizev;
}

double GasGrain::surfaceH2FormationRateCoeff(const GasModule::GrainInterface& gInterface,
                                             double Tgas)
{
	// See Rollig et al. (2013) appendix C + erratum of 2002 Cazaux and Tielens paper
	double total{0};

	for (const auto& pop : *gInterface.populationv())
	{
		const Array& coeffPerGrainPerHPerSizev =
		                surfaceH2FormationRateCoeffPerSize(pop, Tgas);
		total += (pop.densityv() * coeffPerGrainPerHPerSizev).sum();
	}
	return total;
}

Array GasGrain::surfaceH2FormationHeatPerSize(const GasModule::GrainInterface::Population& pop,
                                              double Tgas, double nH)
{
	const Array& coeffPerGrainPerHPerSizev = surfaceH2FormationRateCoeffPerSize(pop, Tgas);
	return coeffPerGrainPerHPerSizev * nH * pop.type()->heatPerH2();
}
