#include "LevelSolution.hpp"
#include "NLevel.hpp"

Array LevelSolution::emissivityv(const Array& eFrequencyv) const
{
	Array total(eFrequencyv.size());
	_levelCoefficients->forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = _levelCoefficients->lineIntensityFactor(upper, lower,
		                                                        _nv(upper));
		LineProfile lp = _levelCoefficients->lineProfile(upper, lower, _t, _cvv);
		lp.addToBinned(eFrequencyv, total, factor);
	});
	return total;
}

Array LevelSolution::opacityv(const Array& oFrequencyv) const
{
	Array total(oFrequencyv.size());
	_levelCoefficients->forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = _levelCoefficients->lineOpacityFactor(upper, lower, _nv(upper),
		                                                      _nv(lower));
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
	for (size_t up = 0; up < _levelCoefficients->numLv(); up++)
	{
		for (size_t lo = 0; lo < up; lo++)
		{
			double cul = _cvv(up, lo);
			double clu = _cvv(lo, up);
			if (cul > 0)
				total += (ev(up) - ev(lo)) * (cul * _nv(up) - clu * _nv(lo));
		}
	}
	return total;
}
