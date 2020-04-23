#ifndef CORE_BIGH2MODEL_HPP
#define CORE_BIGH2MODEL_HPP

#include "H2Data.hpp"
#include "H2Model.hpp"

namespace GasModule
{
    class BigH2Model : public H2Model
    {
    public:
        BigH2Model(const H2Data* h2Data, const Spectrum* meanIntensity);
        void solve(double n, const CollisionParameters& cp, double h2form = 0) override;
        double dissociationRate() const override;
        double dissociationHeating() const override;
        double netHeating() const override;
        double orthoPara() const override;
        Array emissivityv(const Array& eFrequencyv) const override;
        Array opacityv(const Array& oFrequencyv) const override;
        bool hasLevels() const override { return true; }

        /** Read acces to the level solution, should one wish to inspect the solution in more
            detail */
        const LevelSolution* levelSolution() const override { return &_levelSolution; }

        void extraDiagnostics(GasDiagnostics&) const override;

    private:
        /** For each level of X, calculate the direct radiative dissociation rate [s-1] by
            integrating the cross section over the appropriate frequency range. If heatRate is
            true, the integrand is multiplied by (hnu - threshold), so that the heating rate [erg
            s-1] due to this process is given instead. */
        EVector directDissociationIntegralv(bool heatRate = false) const;

        const H2Data* _h2Data;
        const Spectrum* _meanIntensity;
        LevelSolution _levelSolution;
        double _n{0.};

        // values that depend on radiation field only
        EVector _directDissociationRatev;
        EVector _directDissociationHeatv;
    };
}

#endif  // CORE_BIGH2MODEL_HPP
