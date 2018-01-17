#ifndef GASMODULE_GIT_SRC_LINEPROFILE_H
#define GASMODULE_GIT_SRC_LINEPROFILE_H

#include "Array.h"

/** A set of tools to handle line profile in an efficient way. */
class LineProfile
{
public:
	LineProfile(double center, double sigma_gauss, double halfWidht_lorenz);
	// Evaluate the line profile at the given value (frequency) x
	double operator()(double nu) const;
	double center() const { return _center; }

	/** Adds the contribution of a single line to the given spectrum. This way we can stop
	    evaluating the voigt function for the line once the contribution to the total
	    spectrum drops below a chosen threshold. 'factor' is the factor by which the line
	    profile should be multiplied before its values are added to the spectrum. */
	void addToSpectrum(const Array& frequencyv, Array& spectrumv, double factor) const;

	/** Integrate the product of the line profile and the given spectrum in an efficient
	    way. The line profile is only evaluated if the expected contribution of a frequency
	    point to the total exceeds a minimum threshold. */
	double integrateSpectrum(const Array& frequencyv, const Array& spectrumv) const;

private:
	double _center, _sigma_gauss;
	double _one_sqrt2sigma;
	double _a;
};
#endif /* GASMODULE_GIT_SRC_LINEPROFILE_H */
