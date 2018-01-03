#ifndef GASMODULE_GITGASMODULE_GIT_SRC_NLEVEL_H_
#define GASMODULE_GITGASMODULE_GIT_SRC_NLEVEL_H_

#include "Array.h"
#include "EigenAliases.h"
#include "SpecialFunctions.h"

#include <array>
#include <memory>

class LevelDataProvider;

// TODO: use the correct mass for thermal velocity (line width)

/** This class contains a generic implementation for calculating the statistical equilibrium of
    a system with a certain amount of levels (private part), and some properties that result
    from it (public part). NLevel can be used as-is, if this generic (lines only) implementation
    is sufficient, or a subclass of NLevel can be used if extra effects are desired. (e.g. the
    subclass HydrogenLevels provides some extensions to add the two-photon continuum.)

    Example: To simulate the toy two-level system for CII, create an NLevel object taking a
    TwoLevelHardcoded* as constructor argument. To simulate Hydrogen including its two-photon
    continuum, create an object of the HydrogenCalculator subclass, using a HydrogenDataProvider
    class of choice.

    It contains data members to store all the constants needed for these calculations: The
    number of levels, the energies of these levels and their degeneracies, and the spontaneous
    transition rates between them. Since an object of this clas is to be reused many times, or
    by multiple threads at the same time, this class is designed to be stateless: notice how
    almost all of the methods are const.

    However, the calculation of the statisitical equilibrium does need some non-constant
    variables such as the temperature, population densities or the current collision rates.
    These are stored into an object of the type NLevel::Solution, a struct defined in this
    class. A Solution object can be obtained by calling the solveBalance function. By giving
    this object to one of the public functions (e.g. emissivityv), various properties that
    depend on the densities, temperature, collision rates... can be calculated.

    All of the necessary data for the equilibrium calculation is obtained through a class called
    LevelDataProvider. The latter is an abstract class, and the type of atom/molecule that is
    simulated by the implementation in NLevel class will depend on the subclass of
    LevelDataProvider that is used. Each subclass of NLevel will determine what kind of
    LevelDataProvider is used during construction. Some functions of LevelDataProvider are
    called during construction of the NLevel object to set the various constants. In contrast,
    the functions of the LevelDataProvider for the collision coefficients are called at every
    invocation of solveBalance, as these depend on the temperature and collision partner
    densities.

    Subclasses are allowed to override some of the implementations in this class, but they must
    use these same data members. That's why a set of protected getters is provided. */
class NLevel
{
public:
	/** A subclass needs to pass a pointer to a LevelDataProvider object, or the caller can
	    provide one if the generic implementation is sufficient. Of course, it needs to be
	    made sure that this object exists for the lifetime of this class instance. The
	    number of levels, energies, multiplicities, and transition coefficients are filled
	    in immediately using the data provider referred to by this pointer. */
	NLevel(std::shared_ptr<const LevelDataProvider> ldp, const Array& frequencyv);

	virtual ~NLevel();

	/** Getter for the frequency grid the calculations are performed on. All of the input
	    (ouput) spectra must have (will have) the same frequency points. */
	const Array& frequencyv() const { return _frequencyv; }

	/** Ouputs some properties about the different line transitions taken into account by
	    this instance of NLevel. The results for the number of lines, their frequencies and
	    their natural width are returned by reference. */
	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	/** Variables which are changed at every invocation of solveBalance or throughout the
	    calculation, will be passed around using this struct (instead of storing them as
	    members, which is not threadsafe. Note that Solution objects are not interchangeable
	    between NLevel instances, as the number of levels and their properties can be
	    different.

	    TODO: Reconsider how we can prevent breaking abstraction here. Ideally, we'd want a
	    mechanism so that Solutions can only be used with the object that created them.
	    Store pointer to parent? */
	typedef struct Solution
	{
		/* The total density and temperature of the ensemble of atoms/molecules for
		   which the solution was calculated */
		double n, T;

		/* The density of each level population (cm-3) */
		EVector nv;

		/* The induced radiative transition rates and collisional transition rates for
		   this configuration. These are needed to calculate for example the line
		   broadening. */
		EMatrix bpvv, cvv;
	} Solution;

	/** Calculates the level populations for a certain electron temperature and isrf. The
	    temperature, and the densities of the gas species are needed to determine the
	    collisional transition rates and recombination rates. A struct of the type Solution
	    is returned, with all of its members correctly filled in. The collision and
	    recombination data to perform the calculation are obtained from the
	    LevelDataProvider, while the induced transitions rates (B coefficients * line power)
	    are derived from the given specific intensity.

	    External processes can influence the level balance by their contribution to the sink
	    and source terms. Examples are: level-specific formation (think in the chemical
	    network), ionization from specific levels, dissociation from specific levels... Some
	    knowledge about the nature of the levels will be necessary. Subclasses of NLevel can
	    provide an interface to whatever that knowledge may be (usually a mapping from
	    quantum numbers to level index), but the general NLevel implementation will remain
	    completely oblivious. It is up to the client to correctly construct the source and
	    sink terms using this information.

	    More specifically, this function calculates all the matrices for the statistical
	    equilibrium equations, and then calls @c solveRateEquations() do the actual
	    calculation. This function should be generic, while the latter is allowed to have a
	    specialized implementation per subclass. */
	Solution solveBalance(double density, const EVector& speciesNv, double T,
	                      const Array& specificIntensityv, const EVector& sourcev,
	                      const EVector& sinkv) const;

	/** Calculates the level populations using a simple Boltzman LTE equation. */
	Solution solveLTE(double density, const EVector& speciesNv, double T,
	                  const Array& specificIntensityv) const;

	Solution solveZero(double T) const;

	/** The total emitted spectrum by the system of levels. The default implementation gives
	    just the line emission, but subclasses can override it to add extra contributions,
	    such as two-photon continua. */
	virtual Array emissivityv(const Solution& s) const;

	/** The spectrum emitted by the line transitions, expressed as the emission coefficient
	    j_nu f * (erg/cm3/s/Hz). */
	Array lineEmissivityv(const Solution& s) const;

	/** The opacity alpha_nu, equivalent to kappaRho for dust (cm-1). The default
	    implementation gives just the line opacity, but subclasses can override it to add
	    extra contributions, such as absorption cross sections per level. */
	virtual Array opacityv(const Solution& s) const;
	Array lineOpacityv(const Solution& s) const;

	/** Heating rate due to collisional de-excitation (ergs / s / cm3). */
	double heating(const Solution& s) const;

	/** Cooling rate due to collisional excitation. */
	double cooling(const Solution& s) const;

	/** Return the number of levels in the solution */
	size_t numLv() const { return _numLv; }

protected:
	/** Following the notation of the gasPhysics document, construct the rate matrix M_ij =
	    A_ji + B_ji * P_ji + C_ji. Set up F and b using M_ij and the external source term
	    ne*np*alpha_i, due to recombination. Returns the solution as a vector [cm-3].

	    This is a basic implementation which calls a linear solver from the Eigen library. I
	    have make this a virtual function, to make it possible for subclasses to have more
	    specialized algorithm, based on beforehand knowledge about the coefficients. For the
	    H2 model for example, the calculation can be done faster and more precisely by using
	    an iterative approach based on the fact that there is no transition data between and
	    within the electronically excited levels. */
	virtual EVector solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
	                                   const EVector& sourcev, const EVector& sinkv,
	                                   int chooseConsvEq) const;

	EVector solveBoltzmanEquations(double T) const;

private:
	/** Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from
	    the Aij) and Pij is the line power (isrf integrated over the line profile). The
	    units of Bij and Pij are often different in the literature and other codes, but
	    their product should always have units [s-1]. */
	EMatrix prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
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

	/** A wrapper around the voigt function (just some utility). */
	double voigt(double sigma_nu, double halfWidth, double deltaNu) const;

	/** Calculates the Voigt profile for a certain line, using the frequency grid supplied
	    at construction and the temperature and collision rates contained in the Solution
	    struct. */
	Array lineProfile(size_t upper, size_t lower, const Solution& s) const;

	/** Or when the full solution is not yet known, and hence a Solution object is not yet
	    available. */
	Array lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const;

	/** Adds the contribution of a single line to the given spectrum. This way we can stop
	    evaluating the voigt function for the line once the contribution to the total
	    spectrum drops below a chosen threshold. 'factor' is the factor by which the line
	    profile should be multiplied before its values are added to the spectrum. */
	void addLine(Array& spectrumv, size_t upper, size_t lower, double factor, double T,
	             const EMatrix& Cvv) const;

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

	/** Creates the matrix Mij as defined in my notes (basically, Mij = Aji + BPji + Cji). I
	    reused this formula once, so I in a function it goes. */
	EMatrix netTransitionRate(const EMatrix& BPvv, const EMatrix& Cvv) const;

private:
	/** Variables which are the same for all invocations of solveBalance are stored as
	    member. They are set during construction. */

	/* A polymorphic LevelDataProvider. The specific subclass that this data member is
	   initialized with depends on the subclass. */
	std::shared_ptr<const LevelDataProvider> _ldp;

	/* Wavelength grid */
	const Array& _frequencyv;

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

	/** Tabulated voigt function. */
	SpecialFunctions::LookupTable2D _voigt;
};

#endif /* _NLEVEL_H_ */
