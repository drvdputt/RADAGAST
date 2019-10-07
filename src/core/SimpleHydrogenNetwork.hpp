#ifndef CORE_SIMPLEHYDROGENNETWORK_H_
#define CORE_SIMPLEHYDROGENNETWORK_H_

#include "ChemicalNetwork.hpp"

class SimpleHydrogenNetwork : public ChemicalNetwork
{
public:
	/** Sets up a chemical network describing only the reactions between e-, H+, H and
	    H2. */
	SimpleHydrogenNetwork();

	/** Provides the reactions coefficients for each reaction in this network. Specifically, they are, in order:

	    - photoionization of H [s-1]

	    - ionization by collission with electron [cm3 s-1]

	    - radiative recombination [cm3 s-1]

	    - dissociation of H2 after excitation [s-1]

	    - H2 formation (twice the H2 formation rate, as we use 1 H -> 0.5 H2 to make the
              reaction scale linearly with nH) [s-1]
	*/
	EVector rateCoeffv(double T, const Spectrum& specificIntensity,
	                   double kDissFromH2Levels, double kH2FormationGrain) const override;
};

#endif /* CORE_SIMPLEHYDROGENNETWORK_H_ */
