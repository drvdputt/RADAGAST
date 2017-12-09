#ifndef GASMODULE_GIT_SRC_SPECTRUM_H_
#define GASMODULE_GIT_SRC_SPECTRUM_H_

#include "TemplatedUtils.h"

class Spectrum
{
public:
	/** Create an empty object. Can be used as a placeholder in situations where this seems
	    handy. Example: An absorption spectrum for each level of a molecule, but there is no
	    data for some of the levels. */
	Spectrum() : _hasData{false} {};

	/** Create a spectrum, meaning whatever quantity (emissivity, absorption, cross section)
	    as a function of frequency. */
	Spectrum(const Array& frequencyv, const Array& valuev)
	                : _hasData{true}, _frequencyv{frequencyv}, _valuev{valuev} {};

	/** Returns false if no data is available. */
	bool hasData() const { return _hasData; }

	/** Interpolates the data to the given frequency. */
	double evaluate(double frequency) const
	{
		return TemplatedUtils::evaluateLinInterpf(frequency, _frequencyv, _valuev);
	}

private:
	bool _hasData;
	Array _frequencyv, _valuev;
};
#endif /* GASMODULE_GIT_SRC_SPECTRUM_H_ */
