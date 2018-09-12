#ifndef GASMODULE_GIT_SRC_NLEVEL_H_
#define GASMODULE_GIT_SRC_NLEVEL_H_

#include "EigenAliases.h"
#include "LineProfile.h"
#include "SpecialFunctions.h"
#include "Spectrum.h"

#include <array>
#include <memory>

struct GasStruct;
class LevelDataProvider;

/** This class contains a generic implementation for calculating the transitions rates
    determining the statistical equilibrium of a system with a certain amount of levels, and
    some properties that can be derived from the solution. NLevel can be used as-is, if this
    generic (lines only) implementation is sufficient, or a subclass of NLevel can be used if
    extra effects are desired. (e.g. the subclass HydrogenLevels provides some extensions to add
    the two-photon continuum.)

    Examples: 
    
    - To simulate the toy two-level system for CII, create an NLevel object taking a
      TwoLevelHardcoded* as constructor argument.
    
    - To simulate Hydrogen including its two-photon continuum, create an object of the
      HydrogenLevels subclass, using a HydrogenDataProvider class of choice.

    This class stores information about a system of energy levels: the number of levels, the
    energies of these levels and their degeneracies, and the spontaneous transition rates
    between them.

    However, the calculation of the statistical equilibrium and derived properties such as the
    line spectrum does need some non-constant variables such as the temperature, the current
    population densities or the current collision rates. These are stored into an object of the
    type NLevel::Solution, a struct defined in this class. Various properties that depend on the
    densities, temperature, collision rates... can be calculated by passing such an instance
    around. */
class NLevel
{
public:
	/** All of the necessary data is obtained through a class called LevelDataProvider. The
	    latter is an abstract class, and the type of atom/molecule that is simulated by the
	    implementation in NLevel class will depend on the subclass/configuration of the
	    LevelDataProvider instance that is used. Some subclasses of NLevel have extra
	    features, which need more advanced information. This will be made clear in their
	    constructor, either by requiring extra arguments, or by requiring a specific
	    subclass of LevelDataProvider. Constant data is extracted from the LevelDataProvider
	    at construction. In contrast, the functions of the LevelDataProvider for the
	    collision coefficients are called at every calculation of the transition matrix, as
	    these depend on the temperature and collision partner densities.

	    Additionally, the mass of the particle is needed. It determines the thermal velocity
	    of the particle described by this NLevel system, and hence the thermal broadening of
	    the lines. */
	NLevel(std::shared_ptr<const LevelDataProvider> ldp, double mass);

	virtual ~NLevel();

	/** Ouputs some properties about the different line transitions taken into account by
	    this instance of NLevel. The results for the number of lines, their frequencies and
	    their natural width are returned by reference. */
	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	/** Construct the rate matrix T_ij for the given radiation field and gas properties.
	    This is the sum of the spontaneous (Aij), induced (Bij) and collisional (Cij)
	    transitions. The collision data are obtained from the LevelDataProvider, while the
	    induced transitions rates (B coefficients * line power) are derived from the given
	    specific intensity. Optionally, the collision coefficients can be returned
	    separately by pointer, so they don't have to be calculated again later. [s-1] */
	EMatrix totalTransitionRatesvv(const Spectrum& specificIntensity,
	                               const GasStruct& gas, EMatrix* cvv_p = nullptr) const;

	/** Note that Solution objects are not interchangeable between NLevel instances, as the
	    number of levels and their properties can be different. Ideally, we'd want a
	    mechanism so that Solutions can only be used with the object that created them.
	    (Store pointer to parent?). This breaks abstraction a bit. Maybe there are other
	    approaches, but right now we do no checking as that would be tedious. */
	typedef struct Solution
	{
		/* The total density and temperature of the ensemble of atoms/molecules for
		   which the solution was calculated */
		double n, T;

		/* The density of each level population (cm-3) */
		EVector nv;

		/* The collisional transition rates for this configuration. These are needed to
		   calculate for example the line broadening. */
		EMatrix cvv;
	} Solution;

	/** Calculates the level populations using a simple Boltzman LTE equation. */
	Solution solveLTE(double density, const Spectrum& specificIntensity,
	                  const GasStruct& gas) const;

	Solution solveZero(double T) const;

	/** The total emitted spectrum by the system of levels. The default implementation gives
	    just the line emission, but subclasses can override it to add extra contributions,
	    such as two-photon continua. */
	virtual Array emissivityv(const Solution& s, const Array& eFrequencyv) const;

	/** The spectrum emitted by the line transitions, expressed as the emission coefficient
	    j_nu f * (erg/cm3/s/Hz). */
	Array lineEmissivityv(const Solution& s, const Array& eFrequencyv) const;

	/** The opacity alpha_nu, equivalent to kappaRho for dust (cm-1). The default
	    implementation gives just the line opacity, but subclasses can override it to add
	    extra contributions, such as absorption cross sections per level. */
	virtual Array opacityv(const Solution& s, const Array& oFrequencyv) const;
	Array lineOpacityv(const Solution& s, const Array& oFrequencyv) const;

	/** Heating rate due to collisional de-excitation (ergs / s / cm3). */
	double heating(const Solution& s) const;

	/** Cooling rate due to collisional excitation. */
	double cooling(const Solution& s) const;

	/** Return the number of levels in the solution */
	size_t numLv() const { return _numLv; }

	EVector solveBoltzmanEquations(double T) const;

private:
	/** Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from
	    the Aij) and Pij is the line power, i.e. the radiation field integrated over the
	    line profile. The temperature and collsional coefficients are needed, because these
	    influence the shape of the line profile. The units of Bij and Pij are often
	    different in the literature and other codes (it depends on the units used for the
	    radiation field), but their product should always have units [s-1]. */
	EMatrix prepareAbsorptionMatrix(const Spectrum& specificIntensity, double T,
	                                const EMatrix& Cvv) const;

	/** Abstraction of the loop over all lines. Executes thingWithLine for all combinations
	    upper > lower that have _Avv(upper, lower) > 0. If the levels are sorted, and all
	    downward transitions have line activity, then this function will loop over all
	    elements of the lower triangle of the level matrix. */
	void forActiveLinesDo(std::function<void(size_t ini, size_t fin)> thing) const;

	/** Calculates the intensity of a specific line [erg / s / cm2]. Multiplying with the
	    line profile will yield the specific intensity. */
	double lineIntensityFactor(size_t upper, size_t lower, const Solution& s) const;

	/** Computes the opacity of a line, not yet multiplied with the line profile. */
	double lineOpacityFactor(size_t upper, size_t lower, const Solution& s) const;

	/** Return a line profile object that can be used to calculate the (normalized to 1)
	    line profile of the "upper-lower" line. Uses the temperature and collision rates
	    stored in the provided Solution struct. */
	LineProfile lineProfile(size_t upper, size_t lower, const Solution& s) const;

	/** Return a line profile object that can be used to calculate the (normalized to 1)
	    line profile of the "upper-lower" line. The natural line width (the lorenzian
	    contribution) is calculated from the total decay rate due to both spontaneous and
	    collisional transitions contained in _avv and Cvv, respectively. The thermal line
	    width (the gaussian contribution) is calculated from the given temperature and the
	    mass of the particle. */
	LineProfile lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const;

protected:
	/** A number of protected getters are provided, so the subclasses can make use of these
	    coefficients. */
	EVector ev() const { return _ev; }
	double ev(size_t i) const { return _ev(i); }

	EVector gv() const { return _gv; }
	double gv(size_t i) const { return _gv(i); }

	EMatrix avv() const { return _avv; }
	double avv(size_t upper, size_t lower) { return _avv(upper, lower); }

	EMatrix extraAvv() const { return _extraAvv; }
	double extraAvv(size_t upper, size_t lower) const { return _extraAvv(upper, lower); }

private:
	/** Variables which are the same for all invocations of solveBalance are stored as
	    member. They are set during construction. */

	/* A polymorphic LevelDataProvider. The specific subclass that this data member is
	   initialized with depends on the subclass. */
	std::shared_ptr<const LevelDataProvider> _ldp;

	/* Particle mass (important for line width) */
	double _mass;

	/* Energy levels (constant) */
	size_t _numLv{0};
	EVector _ev;

	/* Level degeneracy (constant) */
	EVector _gv;

	/* A matrix (constant, lower triangle, zero diagonal) */
	EMatrix _avv;

	/* Spontaneous transitions that do not produce line photons, but do influence the
	   levels. A prime example is the 2-photon continuum of 2s -> 1s. A2s1 = 2e-6 s-1 for
	   single photon, but is about 8 s-1 for two photons. Maybe this can also be used for
	   some H2 transitions. */
	EMatrix _extraAvv;
};

#endif /* GASMODULE_GIT_SRC_NLEVEL_H_ */
