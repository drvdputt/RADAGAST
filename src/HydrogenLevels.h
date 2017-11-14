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
	HydrogenLevels(std::shared_ptr<const HydrogenDataProvider> hdp, const Array& frequencyv);
	~HydrogenLevels();

	/** This function returns the line emission spectrum + the continuum emitted by the 2s-1s
	    two-photon process. */
	Array emissivityv(const Solution& s) const override;

private:
	/** This function calculates the two-photon continuum using Nussbaumer \& Smutz (1984). */
	Array twoPhotonEmissivityv(const Solution& s) const;

	std::shared_ptr<const HydrogenDataProvider> _hdp;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENLEVELS_H_ */
