#ifndef CORE_BIGH2MODEL_HPP
#define CORE_BIGH2MODEL_HPP

#include "H2Data.hpp"
#include "H2Model.hpp"

class BigH2Model : public H2Model
{
public:
    BigH2Model(const H2Data* h2Data) : _h2Data{h2Data}, _levelSolution(_h2Data) {}
    void solve(double n, const CollisionParameters& cp, const Spectrum& specificIntensity, double h2form = 0) override;
    double dissociationRate(const Spectrum& specificIntensity) const override;
    double dissociationHeating(const Spectrum& specificIntensity) const override;
    double netHeating() const override;
    double orthoPara() const override;
    Array emissivityv(const Array& eFrequencyv) const override;
    Array opacityv(const Array& oFrequencyv) const override;
    bool hasLevels() const override { return true; }

    /** Read acces to the level solution, should one wish to inspect the solution in more
	    detail */
    const LevelSolution* levelSolution() const override { return &_levelSolution; }

private:
    /** Level-resolved dissociation rates. By adding these to the sink terms when solving
	    the statistical equilibrium, the effect of dissociation on the level population can
	    be taken in to account. [s-1] */
    EVector dissociationSinkv(const Spectrum& specificIntensity) const;

    /** For each level of X, calculate the direct radiative dissociation rate [s-1] by
	    integrating the cross section over the appropriate frequency range. If heatRate is
	    true, the integrand is multiplied by (hnu - threshold), so that the heating rate
	    [erg s-1] due to this process is given instead. */
    EVector directDissociationIntegralv(const Spectrum& specificIntensity, bool heatRate = false) const;

    /** Sink term due to the spontaneous dissociation rate. [s-1] */
    const EVector& spontaneousDissociationSinkv() const;

    const H2Data* _h2Data;
    LevelSolution _levelSolution;
    double _n{0.};

    // Keep this allocated, instead of freeing it at the end of solve()
    EMatrix _cvv;
    EMatrix _bpvv;
    EMatrix _tvv;
};

#endif  // CORE_BIGH2MODEL_HPP
