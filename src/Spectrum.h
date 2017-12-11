#ifndef GASMODULE_GIT_SRC_SPECTRUM_H_
#define GASMODULE_GIT_SRC_SPECTRUM_H_

#include "Array.h"
#include "TemplatedUtils.h"

class Spectrum
{
public:
	/** Create an empty object. Can be used as a placeholder in situations where this seems
	    handy. Example: An absorption spectrum for each level of a molecule, but there is no
	    data for some of the levels. */
	Spectrum();

	/** Create a spectrum, meaning whatever quantity (emissivity, absorption, cross section)
	    as a function of frequency. */
	Spectrum(const Array& frequencyv, const Array& valuev);

	/** Returns false if no data is available. */
	bool hasData() const { return _hasData; }

	/** Interpolates the data to the given frequency. */
	double evaluate(double frequency) const;

private:
	bool _hasData;
	double _freqMin{0}, _freqMax{0};
	Array _frequencyv, _valuev;
};
#endif /* GASMODULE_GIT_SRC_SPECTRUM_H_ */
