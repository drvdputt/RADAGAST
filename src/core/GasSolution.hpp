#ifndef CORE_GASSOLUTION_HPP
#define CORE_GASSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "GasState.hpp"
#include "GrainInterface.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"
#include "SpeciesIndex.hpp"

class GasDiagnostics;
class GasInterfaceImpl;
class Spectrum;

/** Objects of this class are created whenever a one of the main functions of GasInterfaceImpl
    are called. This class serves as a workspace, where all the dynamic values (i.e. changing
    while searching for the equilibrium) are stored. Once a GasSolution objects has been
    properly filled in, the public functions can be called to calculate any derived
    quantities.*/
class GasSolution
{
public:
	/** Pointer to gasInterfaceImpl (not sure if really needed, but for now) and inpinging
	    radiation field are constant. Also, need to pass instances of HModel and H2Model
	    here, as they could be abstract (transfer ownership using unique pointer). Working
	    with shared pointer here, while the objects are const. Which means that they can
	    only be changed from outside. Maybe this should be the other way around TODO: think
	    about this. */
	GasSolution(const GasInterfaceImpl* gi, const GasModule::GrainInterface& gri,
	            const Spectrum& specificIntensity,
	            const std::shared_ptr<HModel> hModel,
	            const std::shared_ptr<H2Model> h2Model)
	                : _gasInterfaceImpl{gi}, _grainInterface{gri},
	                  _specificIntensity{specificIntensity}, _hSolution(std::move(hModel)),
	                  _h2Solution(std::move(h2Model))
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
	const GasInterfaceImpl* _gasInterfaceImpl;
	const GasModule::GrainInterface& _grainInterface;
	const Spectrum& _specificIntensity;
	double _t;
	EVector _speciesNv;
	std::shared_ptr<HModel> _hSolution;
	std::shared_ptr<H2Model> _h2Solution;
};

#endif // CORE_GASSOLUTION_HPP
