#ifndef CORE_GASSOLUTION_HPP
#define CORE_GASSOLUTION_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include "H2Model.hpp"
#include "HModel.hpp"

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
	    radiation field are constant */
	GasSolution(const GasInterfaceImpl* gi, const grainInterface& gri,
	            const Spectrum& specificIntensity)
	                : _gasInterfaceImpl{gi}, _grainInterface{gri},
	                  _specificIntensity{specificIntensity}
	{
	}
	/** The temperature */
	double t() const;
	void setT(double t);

	/** The chemistry solution */
	EVector speciesNv() const;
	void setSpeciesNv(EVector nv);
	double nH() const;
	double nH2() const;
	double np() const;

	/** Access to modify H solution  */
	HModel* hSolution();

	/** Access to modify H2 solution */
	H2Model* h2Solution();

	/** The total emissivity per frequency unit, in erg / s / cm^3 / sr / hz */
	Array emisivityv(Array eFrequencyv) const;

	/** The total opacity at each frequency in 1 / cm */
	Array opacityv(Array oFrequencyv) const;

	/** Total cooling */
	double cooling() const;

	/** The total heating, including the grain photoelectric effect, in erg / s / cm^3. */
	double heating() const;

	/** Copies and/or recalculates many diagnostic values, and puts these in the given
	    GasDiagnostics object */
	void fillDiagnostics(GasDiagnostics*) const;

	/** Fills a GasState object with the information contained in this solution */
	void updateGasState(GasState&) const;

private:
	const GasInterfaceImpl* _gasInterfaceImpl;
	const GrainInterface& _grainInterface;
	const Spectrum& _specificIntensity;
	double _t;
	EVector _speciesNv;
	HModel _hSolution;
	H2Model _h2Solution;
};

#endif // CORE_GASSOLUTION_HPP
