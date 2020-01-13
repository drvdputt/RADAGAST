#ifndef CORE_GASSOLUTION_HPP
#define CORE_GASSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "GrainSolution.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"
#include "SpeciesIndex.hpp"

class FreeBound;
class FreeFree;
class GasDiagnostics;
class Spectrum;

/** Objects of this class are created whenever a one of the main functions of GasInterfaceImpl
    are called. This class serves as a workspace, where all the dynamic values (i.e. changing
    while searching for the equilibrium) are stored. There are also functions that update the
    level populations and the grain temperatures and charge distributions. They should typically
    be called after setting a new temperature and a new species density vector. Once the
    equilibria have been calculated, the other public functions can be called to calculate any
    derived quantities. */
class GasSolution
{
public:
    /** Some references to environmental parameters and other models are passed here. For
	    the HModel and H2Model, ownership is transferred (using move) to this object, since
	    their contents change while searching for the solution. When a non-trivial
	    GrainInterface object is passed, this class will keep track of a list of
	    GrainSolution objects. */
    GasSolution(const GasModule::GrainInterface& gri, const Spectrum& specificIntensity,
                const SpeciesIndex* speciesIndex, std::unique_ptr<HModel> hModel, std::unique_ptr<H2Model> h2Model,
                const FreeBound& freeBound, const FreeFree& freeFree);

    GasSolution(GasSolution&&) = default;

    void makeZero();

    /** Solve the level populations for each level model contained here. The formation rate
	    of H2 needs to be passed, because it pumps the H2 level populations. */
    void solveLevels(double formH2 = 0);

    /** Update the charge distribution and temperature of the grains based on the current
	    species densities. Should be called after using setSpeciesNv. */
    void solveGrains();

    /** The radiation field */
    const Spectrum& specificIntensity() const { return _specificIntensity; }

    /** The temperature */
    double t() const { return _t; }

    /** Set a new temperature */
    void setT(double t) { _t = t; }

    /** The chemistry solution */
    const SpeciesVector& speciesVector() const { return _sv; }
    void setSpeciesNv(const EVector& nv) { _sv.setDensities(nv); }
    double nH() const { return _sv.nH(); }
    double nH2() const { return _sv.nH2(); }
    double np() const { return _sv.np(); }
    double ne() const { return _sv.ne(); }

    /** The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
    Array emissivityv(const Array& eFrequencyv) const;

    /** The total opacity at each frequency in 1 / cm */
    Array opacityv(const Array& oFrequencyv) const;

    /** Total cooling */
    double cooling() const;

    /** The total heating, including the grain photoelectric effect, in erg / s / cm^3. */
    double heating() const;

    /** The heating by the grains only (expensive to calculate), minus the cooling by
	    collisions with the grains. Calculated together for efficiency. Optionally returns
	    the individual constributions through the pointer arguments. */
    double grainHeating(double* photoHeat = nullptr, double* collCool = nullptr) const;

    /** Copies and/or recalculates many diagnostic values, and puts these in the given
	    GasDiagnostics object */
    void fillDiagnostics(GasDiagnostics*) const;

    /** Distills the GasSolution object into the necessary information to retrieve opacity
	    and emissivity. Grids on which the opacity and emissivity will be discretized
	    (before being stored in the gas state) need to be provided. */
    GasModule::GasState makeGasState(const Array& oFrequencyv, const Array& eFrequencyv) const;

    /** Get the H2 photodissociation rate from the H2 model [s-1]. */
    double kDissH2Levels() const;

    /** Get the total grain surface H2 formation rate [s-1]. This is per H density unit;
	    multiplying with nH gives [cm-3 s-1]. The contributions of all grain populations
	    (GrainSolution objects) are summed, which use their updated temperature value. */
    double kGrainH2FormationRateCoeff() const;

private:
    std::vector<GrainSolution> _grainSolutionv;
    const Spectrum& _specificIntensity;
    double _t;
    SpeciesVector _sv;
    std::shared_ptr<HModel> _hSolution;
    std::shared_ptr<H2Model> _h2Solution;
    const FreeBound& _freeBound;
    const FreeFree& _freeFree;
};

#endif  // CORE_GASSOLUTION_HPP
