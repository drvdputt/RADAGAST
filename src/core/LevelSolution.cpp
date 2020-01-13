#include "LevelSolution.hpp"
#include "LevelCoefficients.hpp"

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
