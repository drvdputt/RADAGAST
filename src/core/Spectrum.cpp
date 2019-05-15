#include "Spectrum.h"

using namespace std;

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
	if (frequency > _freqMax || frequency < _freqMin || !_hasData)
		return 0;
	else
	{
		double val = TemplatedUtils::evaluateLinInterpf(frequency, _frequencyv,
		                                                _valuev);
		return val;
	}
}

double Spectrum::resolution(double frequency) const
{
	if (!_hasData)
		return 0;

	int index = TemplatedUtils::index(frequency, _frequencyv);
	if (index == 0)
		index++;
	return _frequencyv[index] - _frequencyv[index - 1];
}

double Spectrum::average(double minFreq, double maxFreq) const
{
	if (!_hasData)
		return 0;
	// Integrate over [minFreq, all points inbetween, maxFreq]

	// Integration points
	auto iNuMin = lower_bound(begin(_frequencyv), end(_frequencyv), minFreq);
	auto iNuMax = lower_bound(begin(_frequencyv), end(_frequencyv), maxFreq);
	size_t numFromGrid = distance(iNuMin, iNuMax);

	Array nuv(numFromGrid + 2);
	nuv[0] = minFreq;
	size_t i = 0;
	for (auto nu = iNuMin; nu < iNuMax; nu++)
	{
		nuv[1 + i] = *nu;
		i++;
	}
	nuv[nuv.size() - 1] = maxFreq;

	// Integrand
	auto iValMin = begin(_valuev) + distance(begin(_frequencyv), iNuMin);
	auto iValMax = iValMin + numFromGrid;

	Array integrandv(numFromGrid + 2);
	integrandv[0] = evaluate(minFreq);
	i = 0;
	for (auto val = iValMin; val < iValMax; val++)
	{
		integrandv[1 + i] = *val;
		i++;
	}
	integrandv[integrandv.size() - 1] = evaluate(maxFreq);

	double integral = TemplatedUtils::integrate<double>(nuv, integrandv);
	return integral / (maxFreq - minFreq);
}

Array Spectrum::binned(Array frequencyv) const
{
	Array binnedValuev(frequencyv.size());
	for (size_t i = 0; i < frequencyv.size(); i++)
	{
		// 000 11111 22222 333
		// |--|--.--|--.--|--|
		//i0     1     2     3
		// For the first bin, set the left boundary to freq[0] instead of a halfway point
		double nuMin = i == 0 ? frequencyv[0]
		                      : (frequencyv[i] + frequencyv[i - 1]) / 2.;
		// For the last bin, set the right boundary to freq[n-1] instead of a halfway
		// point
		double nuMax = i == frequencyv.size() - 1
		                               ? frequencyv[frequencyv.size() - 1]
		                               : (frequencyv[i + 1] + frequencyv[i]) / 2.;

		// Take the average over this part of the spectrum
		binnedValuev[i] = average(nuMin, nuMax);
	}
	return binnedValuev;
}

double Spectrum::valMax() const
{
	if (!_hasData)
		return 0;
	auto maxIt = std::max_element(std::begin(_valuev), std::end(_valuev));
	return *maxIt;
}
