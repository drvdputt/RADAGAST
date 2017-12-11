#include "Spectrum.h"

Spectrum::Spectrum() : _hasData{false} {}

Spectrum::Spectrum(const Array& frequencyv, const Array& valuev)
                : _hasData{true}, _frequencyv{frequencyv}, _valuev{valuev}
{
	_freqMin = _frequencyv[0];
	_freqMax = _frequencyv[frequencyv.size() - 1];
}

double Spectrum::evaluate(double frequency) const
{
	if (frequency > _freqMax || frequency < _freqMin)
		return 0;
	else
		return TemplatedUtils::evaluateLinInterpf(frequency, _frequencyv, _valuev);
}
