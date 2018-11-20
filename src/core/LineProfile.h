#ifndef GASMODULE_GIT_SRC_LINEPROFILE_H
#define GASMODULE_GIT_SRC_LINEPROFILE_H

#include "Spectrum.h"

/** A set of tools to handle line profile in an efficient way. */
class LineProfile
{
public:
	LineProfile(double center, double sigma_gauss, double halfWidth_lorentz);
	// Evaluate the line profile at the given value (frequency) x
	double operator()(double nu) const;
	double center() const { return _center; }

	/** Adds the average of the line to the given spectrum in a binned way. The bins over
	    which the averages are taken are defined by the points halfway between the frequency
	    grid points */
	void addToBinned(const Array& frequencyv, Array& binnedSpectrumv, double factor) const;

	/** Adds the contribution of a single line to the given spectrum. This way we can stop
	    evaluating the voigt function for the line once the contribution to the total
	    spectrum drops below a chosen threshold. 'factor' is the factor by which the line
	    profile should be multiplied before its values are added to the spectrum. */
	void addToSpectrum(const Array& frequencyv, Array& spectrumv, double factor) const;

	/** Integrate the product of the line profile and the given spectrum in an efficient
	    way. A recommended set of frequency points for the line (determined by this class)
	    is combined with the frequency points that discretize the given spectrum.

	    The line profile is only evaluated if the expected contribution of a frequency point
	    to the total (using a conservative heuristic) exceeds a minimum threshold.
	    Optionally the maximum of the spectrum can be given, to help with speeding up the
	    calculation (for example when calculation the line integral over the same spectrum
	    for many different lines). This maximum is used to guess when the calculation can be
	    cut off. */
	double integrateSpectrum(const Spectrum& spectrum, double spectrumMax = 0, bool debug=false) const;

private:
	/** Generates a number of points around the line center, spaced according to a power
	    law, within a distance _halfWidth_lorentz * width of the center. */
	Array recommendedFrequencyGrid(int numPoints = 27, double width = 5) const;

	double _center, _sigma_gauss, _halfWidth_lorentz;
	double _one_sqrt2sigma;
	double _a;
};

void test_addToSpectrum();

void test_integrateSpectrum();
#endif /* GASMODULE_GIT_SRC_LINEPROFILE_H */
