#ifndef CORE_LINEPROFILE_H
#define CORE_LINEPROFILE_H

#include "Spectrum.hpp"

namespace RADAGAST
{
    /** A set of tools to handle line profile in an efficient way. */
    class LineProfile
    {
    public:
        LineProfile(double center, double sigma_gauss, double halfWidth_lorentz);

        /** Evaluate the line profile at the given value [Hz-1] */
        double operator()(double nu) const;

        /** Return the center of the line [Hz]*/
        double center() const { return _center; }

        /** Adds the average of the line to the given spectrum in a binned way. The bins over which
            the averages are taken are defined by the points halfway between the frequency grid
            points */
        void addToBinned(const Array& frequencyv, Array& binnedSpectrumv, double factor) const;

        /** Integrate the product of the line profile and the given spectrum in an efficient way. A
            recommended set of frequency points for the line (determined by this class) is combined
            with the frequency points that discretize the given spectrum.

            The line profile is only evaluated if the expected contribution of a frequency point to
            the total (using a conservative heuristic) exceeds a minimum threshold. Optionally the
            maximum of the spectrum can be given, to speed up the calculation (e.g. when
            calculating the line integral over the same spectrum for many different lines). This
            maximum is used to guess when the calculation can be cut off. */
        double integrateSpectrum(const Spectrum& spectrum, double spectrumMax = 0) const;

    private:
        /** Generates a number of points around the line center. For the central points, it is
            attemped keep the vertical spacing more or less constant. For the wings, a linear
            spacing is used. This has been confirmed to work when the line is mostly Gaussian.
            When the line is extremely narrow, compared to the center frequency, only the center
            point is returned, because doing otherwise creates numerical precision problems. */
        Array recommendedFrequencyGrid(int numPoints = 27) const;

        /** Add the line to a single bin of the given discretized spectrum, threating the line
            profile as if it has zero with. The bin chosen is that of which the central
            wavelength is closest to the center of the line. The added value is factor / (bin
            width). Does nothing if the center of the line is outside of the frequency range. */
        void simpleAddToBinned(const Array& frequencyv, Array& binnedSpectrumv, double factor) const;

        double _center, _sigma_gauss, _halfWidth_lorentz;
        double _one_sqrt2sigma;
        double _a;
        // workaround which switches to purely Lorentzian profile when T=0
        bool _lorentzMode{false};
    };
}
#endif  // CORE_LINEPROFILE_H */
