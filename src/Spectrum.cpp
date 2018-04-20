#include "Spectrum.h"

Spectrum::Spectrum() : _hasData{false} {}

Spectrum::Spectrum(const Array& frequencyv, const Array& valuev)
                : _hasData{true}, _frequencyv{frequencyv}, _valuev{valuev}
{
	Error::equalCheck("Size of frequencyv and valuev", _frequencyv.size(), _valuev.size());
	_freqMin = _frequencyv[0];
	_freqMax = _frequencyv[_frequencyv.size() - 1];
}

double Spectrum::evaluate(double frequency) const
{
	if (frequency > _freqMax || frequency < _freqMin)
		return 0;
	else
	{
		double val = TemplatedUtils::evaluateLinInterpf(frequency, _frequencyv,
		                                                _valuev);
		return val;
	}
}

double Spectrum::average(double minFreq, double maxFreq) const
{
	// Find the grid points that lie within the given bounds

	// right of minFreq
	size_t iNuMin = TemplatedUtils::index(minFreq, _frequencyv);
	// right of maxFreq (not included)
	size_t iNuMax = TemplatedUtils::index(maxFreq, _frequencyv);

	// Integration points
	size_t numFromGrid = iNuMax - iNuMin;
	Array nuv(numFromGrid + 2);
	nuv[0] = minFreq;
	for (size_t i = 0; i < numFromGrid; i++)
		nuv[1 + i] = _frequencyv[iNuMin + i];
	nuv[nuv.size() - 1] = maxFreq;

	Array integrandv(nuv.size());
	for (size_t i = 0; i < nuv.size(); i++)
		integrandv[i] = evaluate(nuv[i]);

	double integral = TemplatedUtils::integrate<double>(nuv, integrandv);
	return integral / (maxFreq - minFreq);
}
