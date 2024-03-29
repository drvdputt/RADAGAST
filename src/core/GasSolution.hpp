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

namespace RADAGAST
{
    class FreeBound;
    class FreeFree;
    class GasDiagnostics;
    class Spectrum;

    /** Objects of this class are created whenever a one of the main functions of GasInterfaceImpl
        are called. This class serves as a workspace, where all the dynamic values (i.e. changing
        while searching for the equilibrium) are stored. When each thread has its own GasSolution
        object to work with, the code will be re-entrant.

        The gas temperature T and density vector nv are to be manually set via the available
        setters. The goal of this class is to compute all quantities that can be derived from
        the given (T, nv) state. By having (T, nv) be set manually, this class can also be used
        to recompute more detailed information from just this small number of doubles, which
        might have been discarded by a RADAGAST client for memory efficiency reasons.

        There are functions that update the level populations and the grain temperatures and
        charge distributions. Once the equilibria have been calculated, the other public
        functions can be called to calculate any derived quantities. Some of the latter are not
        const because some caching or calculations using preallocated memory might happen
        somewhere down the composition hierarchy.

        The internal state can be updated consistently by calling the @ updateAll function.

        Concluding, the use case is as follows:

        1. make GasSolution instance

        2. set T and nv

        3. sync all the variables with new T and nv

        4. extract information as desired. The high-level loop of RADAGAST uses this information
           to find a new (T, nv) and goes back to 2.
    */
    class GasSolution
    {
    public:
        /** Some references to environmental parameters and other models are passed here. For the
            HModel and H2Model, ownership is transferred (using move) to this object, since their
            contents change while searching for the solution. When a non-trivial GrainInterface
            object is passed, this class will keep track of a list of GrainSolution objects. */
        GasSolution(const RADAGAST::GrainInterface* gri, const Spectrum* meanIntensity, double fshield,
                    const SpeciesIndex* speciesIndex, std::unique_ptr<HModel> hModel, std::unique_ptr<H2Model> h2Model,
                    const FreeBound* freeBound, const FreeFree* freeFree);

        GasSolution(GasSolution&&) = default;

        void makeZero();

        /** Solve the level populations for each level model contained here. It is recommended
            to call this after solveGrains(), because grains can pump the H2 level
            populations. */
        void solveLevels();

        /** Update the charge distribution and temperature of the grains based on the current
            species densities. Should be called after using setSpeciesNv. */
        void solveGrains();

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

        /** Total cooling, including grain collisions. */
        double cooling() const;

        /** The total heating, including the grain photoelectric effect, in erg / s / cm^3. */
        double heating();

        /** The heating by photoelectric effect on grains. */
        double grainHeating();

        /** The cooling by collisions with grains */
        double grainCooling() const;

        /** Copies and/or recalculates many diagnostic values, and puts these in the given
            GasDiagnostics object */
        void fillDiagnostics(GasDiagnostics*);

        /** Writes the resulting temperature and densities into the given gas state object */
        void setGasState(RADAGAST::GasState&) const;

        /** Get the H2 photodissociation rate from the H2 model [s-1]. */
        double kDissH2Levels() const;

        /** Get the total grain surface H2 formation rate [s-1]. This is per H density unit;
            multiplying with nH gives [cm-3 s-1]. The contributions of all grain populations
            (GrainSolution objects) are summed, which use their updated temperature value. */
        double kGrainH2FormationRateCoeff() const;

        /** Access to the vector of grain solutions (one for each grain population, in the same
            order), for diagnostic purposes. */
        const std::vector<GrainSolution>& grainSolutionv() const { return _grainSolutionv; }

        /** Read access to the H model. Might be replaced later by a set of functions. */
        const HModel* hModel() const { return _hSolution.get(); }

        /** Read access to the H2 model. Might be replaced later by a set of functions. */
        const H2Model* h2Model() const { return _h2Solution.get(); }

    private:
        // externally settable quantities, to be set/updated before solveGrains and solveLevels
        // are called
        double _t;
        SpeciesVector _sv;

        // quantities updated by solveGrains()
        std::vector<GrainSolution> _grainSolutionv;
        double _kGrainH2FormationRateCoeff{0.};

        // workspace for species models
        std::shared_ptr<HModel> _hSolution;
        std::shared_ptr<H2Model> _h2Solution;

        // Radiation field and shielding factor, as supplied by the client code (might be needed
        // in lower-level physics too (probably H2model))
        const Spectrum* _meanIntensity;
        double _fshield{0};

        // Data read in during RADAGAST set up
        const FreeBound* _freeBound;
        const FreeFree* _freeFree;

        // quantities that depend only on the radiation field (constant)
        double _ionHeatPerH;
    };
}
#endif  // CORE_GASSOLUTION_HPP
