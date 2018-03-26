#ifndef GASMODULE_GIT_SRC_H2LEVELS_H_
#define GASMODULE_GIT_SRC_H2LEVELS_H_

#include "NLevel.h"

#include <vector>

struct GasStruct;
class H2FromFiles;

class H2Levels : public NLevel
{
public:
	H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv);
	~H2Levels();

	Array opacityv(const Solution& s) const override;

protected:
	/** An override which is better suited for H2. It scales as a*n^2, where a is the number
	    of iterations, instead of n^3, apparently. I might be interesting to see this with
	    my own eyes. */
	EVector solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
	                           const EVector& sourcev, const EVector& sinkv,
	                           int chooseConsvEq, const GasStruct& gas) const override;

public:
	/** The dissociation rate, both by direct photodissociation and the indirect Solomon
	    process derived from the level population solution. This rate is usable by the
	    chemical network. */
	double dissociationRate(const Solution& s, const Array& specificIntensityv) const;

	/** TODO: There are several processes which have leftover kinetic energy after the
	    dissociation. In a nutshell, these processes convert radiation into kinetic energy.
	    Solomon process dissociation (electronic excitation followed by transition into
	    ground state continuum). Direct radiative dissociation (effect similar to H
	    ionization). I think it should be possible to calculate this from the radiation
	    field, the level populations, and the threshold energy for each level: (photon
	    energy - threshold energy) times the rate. */
	double dissociationHeating(const Solution& s) const;

	/** TODO: Absorption of kinetic energy by collisional dissociation processes. Since this
	    depends on the velocity distribution of the colliding particles, and the energy
	    transferred during the collision, I will need data for this. Have found no
	    candidates yet. */
	double dissociationCooling(const Solution& s) const;

	EVector dissociationSinkv(const Array& specificIntensityv) const;

private:
	/** Sink term due to direct radiative dissociation. Needs radiation field. */
	EVector directDissociationSinkv(const Array& specificIntensityv) const;

	/** TODO: need collisional dissociation contribution. */

	/** Sink term due to the spontaneous dissociation rate. */
	EVector spontaneousDissociationSinkv() const;

private:
	std::shared_ptr<const H2FromFiles> _hff;

	// List of levels that have dissociation cross section data
	std::vector<size_t> _levelsWithCrossSectionv;
};

#endif /* GASMODULE_GIT_SRC_H2LEVELS_H_ */
