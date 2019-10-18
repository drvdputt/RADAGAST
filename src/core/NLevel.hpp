#ifndef CORE_NLEVEL_H_
#define CORE_NLEVEL_H_

#include "EigenAliases.hpp"
#include "LevelSolution.hpp"
#include "LineProfile.hpp"
#include "SpecialFunctions.hpp"
#include "Spectrum.hpp"

#include <array>
#include <memory>

struct GasStruct;
class LevelDataProvider;

/** This class deals with level coefficients. A data class (e.g. H2FromFiles) can inherit from
    this class to re-use the level coefficient infrastructure, if it also wants to read in level
    coefficients. In the typical use case, totalTransitionRatesvv is called, and its output is
    passed to one of the functions in the LevelSolver namespace. 

    The most important function of this class, is prepareAbsorptionMatrix, which integrates over
    each line to calculate the induced transition rates. To facilitate this, a function which
    can generate LineProfile objects is also implemented here, as well as a function to loop
    over all the active (Aij != 0) transitions. */
class LevelCoefficients
{
public:
	/** An atomic (or molecular) mass needs to be passed, as it will influence the thermal
	    broadening of the lines. */
	LevelCoefficients(double mass);
	virtual ~LevelCoefficients();

protected:
	/** For use by subclass during construction (workaround would be a virtual setup()
	    function.) */
	void setConstants(const EVector& ev, const EVector& gv, const EMatrix& avv,
	                  const EMatrix& extraAvv);

public:
	/** Energy of the levels */
	EVector ev() const { return _ev; }

	/** Multiplicity of the levels */
	EVector gv() const { return _gv; }

	/** Spontaneous radiative transition rates between the levels */
	EMatrix avv() const { return _avv; }

	/** Collisional transition rates, calculated by a subclass */
	virtual EMatrix cvv(const GasStruct& gas) const = 0;

	/** Spontaneous transition rates which do not produce line emission (really only used
	    for two-photon continuum). */
	EMatrix extraAvv() const { return _extraAvv; }

	/** Ouputs some properties about the different line transitions. The results for the
	    number of lines, their frequencies [s-1] and their natural widths (decay rate [s-1]
	    / 4 pi) are returned by reference. */
	void lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const;

	/** Construct the rate matrix T_ij for the given radiation field and gas properties.
	    This is the sum of the spontaneous (Aij), induced (Bij) and collisional (Cij)
	    transitions. The collision data are obtained from the LevelDataProvider, while the
	    induced transitions rates (B coefficients * line power) are derived from the given
	    specific intensity. Optionally, the collision coefficients can be returned
	    separately by pointer, so they don't have to be calculated again later (the result
	    is exactly the same as calling cvv() but it is more efficient to calculate the total
	    matrix and cvv at the same time). [s-1] */
	EMatrix totalTransitionRatesvv(const Spectrum& specificIntensity, const GasStruct& gas,
	                               EMatrix* cvv_p = nullptr) const;

	/** Calculates the level populations using a simple Boltzman LTE equation. Also serves
	    as an example of how to properly set up a LevelSolution object. */
	LevelSolution solveLTE(double density, const GasStruct& gas) const;

	LevelSolution solveZero(double T) const;

	/** Return the number of levels in the solution */
	size_t numLv() const { return _numLv; }

	/** The boltzman fractions for the levels, based purely on their energies and the
	    temperatures */
	EVector solveBoltzmanEquations(double T) const;

private:
	/** Create the matrix [Bij*Pij], where Bij are the Einstein B coefficients (derived from
	    the Aij) and Pij is the line power, i.e. the radiation field integrated over the
	    line profile. The temperature and collision coefficients are needed, because these
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

	/** Calculates the integrated emission coefficient of a specific line [erg s-1 cm-3
	    sr-1]. Multiplying with the line profile [Hz-1] will yield the specific intensity
	    [erg s-1 cm-3 sr-1 Hz-1]. */
	double lineIntensityFactor(size_t upper, size_t lower, double nu, double nl) const;

	/** Computes the integrated opacity of a line [cm-1 Hz]. Multiplying with the line
	    profile [Hz-1] will yield the opacity [cm-1] at each frequency. */
	double lineOpacityFactor(size_t upper, size_t lower, double nu, double nl) const;

	/** Return a line profile object that can be used to calculate the (normalized to 1)
	    line profile of the "upper-lower" line. The natural line width (the lorenzian
	    contribution) is calculated from the total decay rate due to both spontaneous and
	    collisional transitions contained in _avv and Cvv, respectively. The thermal line
	    width (the gaussian contribution) is calculated from the given temperature and the
	    mass of the particle. */
	LineProfile lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const;

private:
	double _mass;
	size_t _numLv{0};
	EVector _ev;
	EVector _gv;
	EMatrix _avv;
	EMatrix _extraAvv;
};

#endif /* CORE_NLEVEL_H_ */
