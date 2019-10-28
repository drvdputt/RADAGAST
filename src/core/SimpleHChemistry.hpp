#ifndef CORE_SIMPLEHCHEMISTRY_HPP
#define CORE_SIMPLEHCHEMISTRY_HPP

#include "Chemistry.hpp"

class Spectrum;

class SimpleHChemistry : public Chemistry
{
public:
	SimpleHChemistry();

	/** Calculate the rate coefficients for each reaction. Multiplying with the right
	    density factors of the reactants will give a total rate [cm-3 s-1]. The unit of
	    these coefficients is hence [cm(-3+n) s-1] where n is usually the number of reactant
	    particles. Since the calculation of the H2 formation and destruction rate are
	    complicated, the are expected to be calculated somewhere else, and passed as
	    arguments here. This function simply fills them in into the right spot of the
	    k-vector. */
	EVector rateCoeffv(double T, const Spectrum& specificIntensity,
	                   double kDissFromH2Levels, double kH2FormationGrain) const;
};

#endif // CORE_SIMPLEHCHEMISTRY_HPP
