#ifndef CORE_GASSOLUTION_HPP
#define CORE_GASSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"
#include "SpeciesIndex.hpp"

class FreeBound;
class FreeFree;
class GasDiagnostics;
class Spectrum;

/** Objects of this class are created whenever a one of the main functions of GasInterfaceImpl
    are called. This class serves as a workspace, where all the dynamic values (i.e. changing
    while searching for the equilibrium) are stored. Once a GasSolution objects has been
    properly filled in, the public functions can be called to calculate any derived
    quantities.*/
class GasSolution
{
public:
	/** Some references to environmental parameters and other models are passed here. For
	    the HModel and H2Model, ownership is transferred (using move) to this object, since
	    their contents change while searching for the solution. */
	GasSolution(const GasModule::GrainInterface& gri, const Spectrum& specificIntensity,
	            std::unique_ptr<HModel> hModel, std::unique_ptr<H2Model> h2Model,
	            const FreeBound& freeBound, const FreeFree& freeFree)
	                : _grainInterface{gri}, _specificIntensity{specificIntensity},
	                  _hSolution(std::move(hModel)), _h2Solution(std::move(h2Model)),
	                  _freeBound{freeBound}, _freeFree{freeFree}
	{
	}

	GasSolution(GasSolution&&) = default;

	void makeZero();

	/** Solve the level populations for each level model contained here. The formation rate
	    of H2 needs to be passed, because it pumps the H2 level populations. */
	void solveLevels(double formH2 = 0);

	/** The radiation field */
	const Spectrum& specificIntensity() const { return _specificIntensity; }

	/** The temperature */
	double t() const { return _t; }
	void setT(double t) { _t = t; }

	/** The chemistry solution */
	EVector speciesNv() const { return _speciesNv; }
	void setSpeciesNv(const EVector& nv) { _speciesNv = nv; }
	double nH() const { return _speciesNv(SpeciesIndex::inH()); }
	double nH2() const { return _speciesNv(SpeciesIndex::inH2()); }
	double np() const { return _speciesNv(SpeciesIndex::inp()); }
	double ne() const { return _speciesNv(SpeciesIndex::ine()); }

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
	GasModule::GasState makeGasState(const Array& oFrequencyv,
	                                 const Array& eFrequencyv) const;

	double kDissH2Levels() const;

private:
	const GasModule::GrainInterface& _grainInterface;
	const Spectrum& _specificIntensity;
	double _t;
	EVector _speciesNv;
	std::shared_ptr<HModel> _hSolution;
	std::shared_ptr<H2Model> _h2Solution;
	const FreeBound& _freeBound;
	const FreeFree& _freeFree;
};

#endif // CORE_GASSOLUTION_HPP
