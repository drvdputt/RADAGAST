#include "LevelSolution.hpp"
#include "CollisionParameters.hpp"
#include "LevelCoefficients.hpp"

namespace GasModule
{
    void LevelSolution::updateRates(const Spectrum& meanIntensity, const CollisionParameters& cp)
    {
        _t = cp._t;
        _cvv = _levelCoefficients->cvv(cp);
        _bpvv = _levelCoefficients->prepareAbsorptionMatrix(meanIntensity, _t, _cvv);
    }

    void LevelSolution::updateRates(const CollisionParameters& cp)
    {
        _t = cp._t;
        _cvv = _levelCoefficients->cvv(cp);
        _bpvv = EMatrix::Zero(_levelCoefficients->numLv(), _levelCoefficients->numLv());
    }

    void LevelSolution::setToZero(double T)
    {
        _t = T;
        int numLv = _levelCoefficients->numLv();
        _cvv = EMatrix::Zero(numLv, numLv);
        _bpvv = EMatrix::Zero(numLv, numLv);
        _nv = EVector::Zero(numLv);
    }

    EMatrix LevelSolution::Tvv() const
    {
        return _levelCoefficients->avv() + _levelCoefficients->extraAvv() + _bpvv + _cvv;
    }

    bool LevelSolution::hasBadNv() const
    {
        return _nv.size() == 0 || (_nv.array() < 0).any() || (_nv.array() == 0).all();
    }

    EVector LevelSolution::fv() const
    {
        double sum = _nv.sum();
        if (sum > 0)
            return _nv / sum;
        else
            // return something to helps with e.g. H2 dissociation rate (needs to be non-zero
            // even if H2 density is zero)
            return _levelCoefficients->solveBoltzmanEquations(_t);
    }

    Array LevelSolution::emissivityv(const Array& eFrequencyv) const
    {
        Array total(eFrequencyv.size());
        _levelCoefficients->forActiveLinesDo([&](size_t upper, size_t lower) {
            double factor = _levelCoefficients->lineIntensityFactor(upper, lower, _nv(upper));
            LineProfile lp = _levelCoefficients->lineProfile(upper, lower, _t, _cvv);
            lp.addToBinned(eFrequencyv, total, factor);
        });
        return total;
    }

    Array LevelSolution::opacityv(const Array& oFrequencyv) const
    {
        Array total(oFrequencyv.size());
        _levelCoefficients->forActiveLinesDo([&](size_t upper, size_t lower) {
            double factor = _levelCoefficients->lineOpacityFactor(upper, lower, _nv(upper), _nv(lower));
            LineProfile lp = _levelCoefficients->lineProfile(upper, lower, _t, _cvv);
            // lp.addToSpectrum(oFrequencyv, total, factor);
            lp.addToBinned(oFrequencyv, total, factor);
        });
        return total;
    }

    double LevelSolution::netHeating() const
    {
        double total = 0;
        const EVector& ev = _levelCoefficients->ev();
        for (size_t j = 0; j < _levelCoefficients->numLv(); j++)
        {
            for (size_t i = 0; i < _levelCoefficients->numLv(); i++)
            {
                double Cij = _cvv(i, j);
                if (Cij > 0) total += (ev(i) - ev(j)) * Cij * _nv(i);
            }
        }
        return total;
    }
}
