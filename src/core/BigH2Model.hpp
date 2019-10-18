#ifndef CORE_BIGH2MODEL_HPP
#define CORE_BIGH2MODEL_HPP
#include "H2FromFiles.hpp"
#include "H2Model.hpp"
#include "LevelSolution.hpp"

class BigH2Model : public H2Model
{
public:
	BigH2Model(const H2Data* h2Data) : _h2Data{h2Data}, _levelSolution(_h2Data) {}
	void solve(double n, const GasStruct& gas, const Spectrum& specificIntensity,
	           double h2form = 0) override;
	double dissociationRate(const Spectrum& specificIntensity) const override;
	double orthoPara() const override;
	Array emissivityv(const Array& eFrequencyv) const override;
	Array opacityv(const Array& oFrequencyv) const override;

private:
	/** Level-resolved dissociation rates. By adding these to the sink terms when solving
	    the statistical equilibrium, the effect of dissociation on the level population can
	    be taken in to account. [s-1] */
	EVector dissociationSinkv(const Spectrum& specificIntensity) const;

	/** For each level of X, calculate the direct radiative dissociation rate [s-1] by
	    integrating the cross section over the appropriate frequency range. If heatRate is
	    true, the integrand is multiplied by (hnu - threshold), so that the heating rate
	    [erg s-1] due to this process is given instead. */
	EVector directDissociationIntegralv(const Spectrum& specificIntensityv,
	                                    bool heatRate = false) const;

	/** Sink term due to the spontaneous dissociation rate. [s-1] */
	EVector spontaneousDissociationSinkv() const;

	/** TODO: need collisional dissociation contribution. */

	const H2Data* _h2Data;
	LevelSolution _levelSolution;
	double _n{0.};
};

#endif // CORE_BIGH2MODEL_HPP
