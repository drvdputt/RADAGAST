#include "HModel.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Ionization.hpp"
#include "LevelSolver.hpp"

void HModel::solve(double n, const CollisionParameters& cp, const Spectrum& specificIntensity)
{
    _levelSolution.setT(cp._t);

    if (n <= 0)
    {
        int numLv = _hData->numLv();
        _levelSolution.setCvv(EMatrix::Zero(numLv, numLv));
        _levelSolution.setNv(EVector::Zero(numLv));
        return;
    }

    EMatrix Cvv;
    EMatrix Tvv = _hData->totalTransitionRatesvv(specificIntensity, cp, &Cvv);
    _levelSolution.setCvv(Cvv);

    EVector the_sourcev = sourcev(cp);
    EVector the_sinkv = sinkv(cp);
    DEBUG("Solving levels nH = " << n << std::endl);
    EVector newNv = LevelSolver::statisticalEquilibrium(n, Tvv, the_sourcev, the_sinkv);
    _levelSolution.setNv(newNv);
}

Array HModel::emissivityv(const Array eFrequencyv) const
{
    return _levelSolution.emissivityv(eFrequencyv) + twoPhotonEmissivityv(eFrequencyv);
}

Array HModel::opacityv(const Array oFrequencyv) const
{
    return _levelSolution.opacityv(oFrequencyv);
}

double HModel::netHeating() const
{
    return _levelSolution.netHeating();
}

Array HModel::twoPhotonEmissivityv(const Array& eFrequencyv) const
{
    std::array<size_t, 2> upper_lower = _hData->twoPhotonIndices();
    // This index can mean either the resolved level 2s, or the collapsed level 2
    size_t index2sOr2 = upper_lower[0];
    size_t index1s = upper_lower[1];

    const EVector& gv = _hData->gv();
    const EVector& ev = _hData->ev();

    // The population of the 2s level needs to be guessed when n=2 is collapsed. We can
    // check this by looking at the multiplicity of the level. Since 2s and 1s should both
    // have the the same multiplicity, we can do:
    bool collapsed = gv(index2sOr2) != gv(index1s);

    // If collapsed, assume the population of the 2s level to be 1/4 of the total n=2
    // population.
    double n2s = _levelSolution.nv()(index2sOr2);
    if (collapsed) n2s /= 4.;

    // 1984-Nussbaumer
    // constant factor in eq 3
    double constFactor = Constant::PLANCK / Constant::FPI * n2s;
    double nu0 = (ev(index2sOr2) - ev(index1s)) / Constant::PLANCK;

    // Parameters for eq 2
    const double C = 202.0;  // s-1
    const double alpha = .88;
    const double beta = 1.53;
    const double gam = .8;

    Array result(eFrequencyv.size());
    for (size_t iFreq = 0; eFrequencyv[iFreq] < nu0 && iFreq < eFrequencyv.size(); iFreq++)
    {
        double y = eFrequencyv[iFreq] / nu0;
        double y1miny = y * (1 - y);
        double pow4y1miny_gam = pow(4 * y1miny, gam);
        double Py = C * (y1miny * (1 - pow4y1miny_gam) + alpha * pow(y1miny, beta) * pow4y1miny_gam);
        result[iFreq] = constFactor * y * Py;
    }
    return result;
}

EVector HModel::sourcev(const CollisionParameters& cp) const
{
    EVector result = _hData->recombinationRatev(cp._t);
    double ne = cp._sv.ne();
    double np = cp._sv.np();
    return result * ne * np;
}

EVector HModel::sinkv(const CollisionParameters& cp) const
{
    // TODO: ideally, this calculates the ionization rate from each level, using individual
    // ionization cross sections.

    /* The ionization rate calculation makes no distinction between the levels.  When
	   the upper level population is small, and its decay rate is large, the second term
	   doesn't really matter. Therefore, we choose the sink to be the same for each
	   level.  Moreover, total source = total sink so we want sink*n0 + sink*n1 = source
	   => sink = totalsource / n because n0/n + n1/n = 1. */
    double ne = cp._sv.ne();
    double np = cp._sv.np();
    double nH = cp._sv.nH();
    double totalSource = ne * np * Ionization::recombinationRateCoeff(cp._t);
    double sink = totalSource / nH;  // Sink rate per (atom per cm3)
    return EVector::Constant(_hData->numLv(), sink);
}
