#ifndef CORE_SPECTRUM_HPP
#define CORE_SPECTRUM_HPP

#include "Array.hpp"
#include "TemplatedUtils.hpp"

namespace RADAGAST
{
    class Spectrum
    {
    public:
        /** Create an empty object. Can be used as a placeholder in situations where this seems
            handy. Example: An absorption spectrum for each level of a molecule, but there is no
            data for some of the levels. */
        Spectrum();

        /** Create a spectrum, meaning whatever quantity (emissivity, absorption, cross section) as
            a function of frequency. */
        Spectrum(const Array& frequencyv, const Array& valuev);

        /** Returns false if no data is available. */
        bool hasData() const { return _hasData; }

        /** Returns the number of frequency points */
        int numPoints() const { return _frequencyv.size(); }

        /** Interpolates the data to the given frequency. */
        double evaluate(double frequency) const;

        /** Distance between data points between which the given frequency falls. If the given
            frequency is outside of the frequency grid, then the distance between the first or last
            pair of points in the frequency grid is given. Returns 0 if the spectrum is empty. */
        double resolution(double frequency) const;

        /** Calculate the average of the spectrum over a given interval. The points used to perform
            the integration are the points on which the spectrum is discretized, with minFreq and
            maxFreq added extra to fix the range. */
        double average(double minFreq, double maxFreq) const;

        /** Bin this spectrum onto another frequency grid. The argument given is a list of
            frequencies, on which this spectrum should be discretized. Every value of the resulting
            Array, is this spectrum averaged over a bin surrounding each frequency point given in
            frequencyv. The bin edges are put halfway between the frequency points. The first and
            last bin are special cases: they go from frequencyv[0] to the first intermediate point,
            and from the last intermediate point to frequencyv[n-1]. */
        Array binned(Array frequencyv) const;

        double freqMax() const { return _freqMax; }
        double freqMin() const { return _freqMin; }
        double valMax() const;
        double valMin() const;

        const Array& frequencyv() const { return _frequencyv; }
        const Array& valuev() const { return _valuev; }

    private:
        bool _hasData;
        double _freqMin{0};
        double _freqMax{0};
        Array _frequencyv, _valuev;
    };
}
#endif  // CORE_SPECTRUM_HPP
