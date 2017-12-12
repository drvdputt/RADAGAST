#ifndef GASMODULE_GIT_SRC_HYDROGENLEVELS_H_
#define GASMODULE_GIT_SRC_HYDROGENLEVELS_H_

#include "NLevel.h"

class HydrogenDataProvider;

/** This class expands \c NLevel with the functionality necessary to simulate a Hydrogen atom. In
    practice, this includes:

    - Providing a constructor which takes \c HydrogenDataProvider as an argument, instead of a
      regular \c LevelDataProvider.

    - An override for \c emissivityv, which adds the two-photon continuum to the standard
      implementation which only returns the line emission. */
class HydrogenLevels : public NLevel
{
public:
	/** Work with a shared pointer here, so that a pointer can be given both to the derived
	    class and the base class. */
	HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp,
	               const Array& frequencyv);
	~HydrogenLevels();

	/** This function returns the line emission spectrum + the continuum emitted by the 2s-1s
	    two-photon process. */
	Array emissivityv(const Solution& s) const override;

	/** Returns a vector containing the source terms for the equilibrium equations, such as the
	    partial recombination rates into each level. Note that this function is separated from
	    the ionization balance calculation, as there only the total recombination rate matters.
	    In this case, this function returns the partial recombination rate into each level,
	    interpolated from some formulae found in 2015-Raga (for levels 3-5) , and Draine's book
	    (for levels 1, 2). The l-resolved recombination rates are just weighed by the degeneracy
	    of the level, for levels 3, 4 and 5. TODO: Use better data here. [cm-3 s-1] */
	EVector sourcev(double T, const EVector& speciesNv) const;

	/** Produces the sink term to be used by the equilibrium equations. In this case, hydrogen
	    disappears from the level populations because it's being ionized. In the current
	    implementation, all the ionization is assumed to be drawn equally from all levels. TODO:
	    add the effects of H2 formation in here?. Take care of this using actual ionization
	    cross sections? [s-1] */
	EVector sinkv(double T, double n, const EVector& speciesNv) const;

private:
	/** This function calculates the two-photon continuum using Nussbaumer \& Smutz (1984). */
	Array twoPhotonEmissivityv(const Solution& s) const;

	std::shared_ptr<const HydrogenDataProvider> _hdp;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENLEVELS_H_ */
