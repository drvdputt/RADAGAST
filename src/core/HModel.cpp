#include "HModel.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "EigenAliases.hpp"
#include "Ionization.hpp"
#include "LevelSolver.hpp"
#include "Spectrum.hpp"
#include "TwoPhoton.hpp"

namespace GasModule
{
    HModel::HModel(const HData* hData, const Spectrum* meanIntensity)
        : _hData{hData}, _meanIntensity{meanIntensity}, _levelSolution(_hData)
    {
        _ionizationRate = Ionization::photoRateCoeff(*meanIntensity);
    }

    void HModel::solve(double n, const CollisionParameters& cp)
    {
        if (n <= 0)
        {
            _levelSolution.setToZero(cp._t);
            return;
        }
        DEBUG("Solving levels nH = " << n << std::endl);

        _levelSolution.updateRates(*_meanIntensity, cp);
        EVector newNv = LevelSolver::statisticalEquilibrium(n, _levelSolution.Tvv(), sourcev(cp), sinkv());
        _levelSolution.setNv(newNv);
    }

    Array HModel::emissivityv(const Array eFrequencyv) const
    {
        return _levelSolution.emissivityv(eFrequencyv) + twoPhotonEmissivityv(eFrequencyv);
    }

    Array HModel::opacityv(const Array oFrequencyv) const { return _levelSolution.opacityv(oFrequencyv); }

    double HModel::netHeating() const { return _levelSolution.netHeating(); }

    double HModel::n2s() const
    {
        auto upper_lower = _hData->twoPhotonIndices();
        // This index can mean either the resolved level 2s, or the collapsed level 2
        size_t index2sOr2 = upper_lower[0];
        size_t index1s = upper_lower[1];

        // The population of the 2s level needs to be guessed when n=2 is collapsed. We can check
        // this by looking at the multiplicity of the level. Since 2s and 1s should both have the
        // the same multiplicity, we can do:
        const EVector& gv = _hData->gv();
        bool collapsed = gv(index2sOr2) != gv(index1s);

        // If collapsed, assume the population of the 2s level to be 1/4 of the total n=2
        // population.
        if (collapsed)
            return 0.25 * _levelSolution.nv()(index2sOr2);
        else
            return _levelSolution.nv()(index2sOr2);
    }

    Array HModel::twoPhotonEmissivityv(const Array& eFrequencyv) const
    {
        return TwoPhoton::emissivityv(eFrequencyv, n2s());
    }

    EVector HModel::sourcev(const CollisionParameters& cp) const
    {
        EVector result = _hData->recombinationRatev(cp._t);
        double ne = cp._sv.ne();
        double np = cp._sv.np();
        return result * ne * np;
    }

    EVector HModel::sinkv() const
    {
        // The ionization rate calculation makes no distinction between the levels. When the
        // upper level population is small, and its decay rate is large, the second term doesn't
        // really matter. Therefore, only add a sink term for the ground state (however, I think
        // the level solver replaces this equation by the conservation equation anyway, so it
        // probably doesn't matter).
        EVector result = EVector::Zero(_hData->numLv());
        result[_hData->index(1, 0)] = _ionizationRate;
        return result;
    }
}
