#include "LineProfile.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Functions.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "TemplatedUtils.hpp"
#include <cmath>
#include <iomanip>
#include <iterator>

using namespace std;

namespace RADAGAST
{
    LineProfile::LineProfile(double center, double sigma_gauss, double halfWidth_lorentz)
        : _center{center}, _sigma_gauss{sigma_gauss}, _halfWidth_lorentz(halfWidth_lorentz)

    {
        _one_sqrt2sigma = M_SQRT1_2 / _sigma_gauss;
        _a = _halfWidth_lorentz * _one_sqrt2sigma;
        if (_sigma_gauss <= 0.) _lorentzMode = true;
    }

    double LineProfile::operator()(double nu) const
    {
        // double x = (nu - _center) * _one_sqrt2sigma;
        // Note that the normalization factor is 1 / sqrt(2 pi) sigma
        // return Functions::voigt(_a, x) / Constant::SQRT2PI / _sigma_gauss;
        // return Functions::pseudoVoigt(nu - _center, _sigma_gauss, _halfWidth_lorentz);

        // use gauss, except for zero gaussian width (to avoid bugs)
        if (_lorentzMode)
            return Functions::lorentz(nu - _center, _halfWidth_lorentz);
        else
            return Functions::gauss(nu - _center, _sigma_gauss);
    }

    Array LineProfile::recommendedFrequencyGrid(int numPoints) const
    {
        // Odd numbers of points make this easier
        if (!(numPoints % 2)) numPoints++;

        // Because then this is the center
        int iCenter = numPoints / 2;

        double fwhmGauss = 2.35 * _sigma_gauss;
        double fwhmLorentz = 2 * _halfWidth_lorentz;
        // Approximation from wikipedia (actual source was is 404)
        double fwhmVoigt = 0.5346 * fwhmLorentz + std::sqrt(0.2166 * fwhmLorentz * fwhmLorentz + fwhmGauss * fwhmGauss);

        double yMax = (*this)(_center);
        // With 2.33 sigma, we get about 0.99 of the gaussian
        double xMax = std::max(2.8 * _sigma_gauss, Functions::lorentz_percentile(0.995, _halfWidth_lorentz));

        // avoid numerical problems for very narrow lines by returning only 1 point
        if (xMax < 1e-15 * _center) return Array(_center, 1);

        double yMin = (*this)(_center + xMax);
        double linearStep = (yMax - yMin) / (iCenter + 1);  // +1 to avoid stepping below 0

        // Since yMax, and hence (y - linearStep) can be higher than the range of the inverse
        // function we use below (i.e. the peak of Lorentz or Gauss), it can run into trouble.
        // Scale y down (equivalent to scaling the gauss up) using this factor.
        double yMaxAllowed;
        if (_lorentzMode)
            yMaxAllowed = Functions::lorentz(0, _halfWidth_lorentz);
        else
            yMaxAllowed = Functions::gauss(0, _sigma_gauss);
        // make sure that maximum of this helper function corresponds to actual maximum of line profile
        double yScale = yMaxAllowed / yMax;

        // Fill in the values
        Array freqv(numPoints);
        freqv[iCenter] = _center;

        double y = yMax;
        // double yStepFactor = pow(untilFactor, 1. / iCenter);f
        double x = 0;
        double previousx = 0;
        // Example: iCenter = 2
        // frequency index:
        // 0 1 2 3 4
        // i:
        // 2 1 0 1 2
        // So i must include iCenter --> use <=
        int i = 0;
        // Within the fwhm, use this fancy algorithm
        for (i = 1; i <= iCenter / 2; i++)
        {
            previousx = x;

            // We want to go down a fixed step in the vertical direction.
            double offset = (y - linearStep) * yScale;
            if (_lorentzMode)
                x = Functions::inverse_lorentz(offset, _halfWidth_lorentz);
            else
                x = Functions::inverse_gauss(offset, _sigma_gauss);

            // If drop was too big, scale down
            double nexty = (*this)(_center + x);
            while (y - nexty > linearStep)
            {
                x -= (x - previousx) / 4;
                nexty = (*this)(_center + x);
            }
            // Do the next step starting from a new y position
            nexty = y - linearStep;
            y = nexty;

            // This is safe because we forced the number of points to be odd
            freqv[iCenter + i] = _center + x;
            freqv[iCenter - i] = _center - x;

            if (x > fwhmVoigt / 2) break;
        }
        // Outside the FWHM, or when having used half of the points, continue linearly, to make
        // sure we get to the end
        int left = iCenter - i;
        double step = (xMax - x) / left;
        for (; i <= iCenter; i++)
        {
            x += step;
            freqv[iCenter + i] = _center + x;
            freqv[iCenter - i] = _center - x;
        }
        // std::cout << "added " << left << " linear points" << '\n';
        return freqv;
    }

    void LineProfile::addToBinned(const Array& frequencyv, Array& binnedSpectrumv, double factor) const
    {
        if (frequencyv.size() < 2) Error::runtime("frequencyv should have at least two values");

        double gridMin = frequencyv[0];
        double gridMax = frequencyv[frequencyv.size() - 1];

        const Array& lineGrid = recommendedFrequencyGrid();

        if (lineGrid.size() == 1)
        {
            // if the line is too narrow to form a proper grid, just add total integrated value
            // (1 * factor / deltanu) to the right bin
            simpleAddToBinned(frequencyv, binnedSpectrumv, factor);
            return;
        }

        // else, use the generated frequency grid
        auto lineBegin = begin(lineGrid);
        auto lineEnd = end(lineGrid);

        // If line is completely outside grid, do nothing
        if (*lineBegin > gridMax || *(lineEnd - 1) < gridMin) return;

        // line wider than grid --> ignore line points not within grid
        while (*lineBegin < gridMin) lineBegin++;
        while (*(lineEnd - 1) > gridMax) lineEnd--;
        double lineMin = *lineBegin;
        double lineMax = *(lineEnd - 1);
        auto numLinePoints = distance(lineBegin, lineEnd);

        // line narrower than grid --> ignore grid points not within span of line
        // (the above code guarantees that lineMax is smaller than the max of frequencyv)
        auto gridBegin = upper_bound(begin(frequencyv), end(frequencyv), lineMin) - 1;
        auto gridEnd = upper_bound(begin(frequencyv), end(frequencyv), lineMax);
        auto numGridPoints = distance(gridBegin, gridEnd);

        // Now merge the bin centers and the line points, which will serve as an
        // integration grid
        Array integrationGridv(numGridPoints + numLinePoints);
        merge(gridBegin, gridEnd, lineBegin, lineEnd, begin(integrationGridv));

        // Calculate the integrand on this grid
        Array integrandv(integrationGridv.size());
        for (size_t i = 0; i < integrandv.size(); i++) integrandv[i] = (*this)(integrationGridv[i]);

        // Go over the bins, and integrate the parts of the line that fall within each bin.
        // Starting with offset (which is were our grid for the integration starts)
        auto beginIndex = distance(begin(frequencyv), gridBegin);
        auto endIndex = distance(begin(frequencyv), gridEnd);
        for (int bin = beginIndex; bin < endIndex; bin++)
        {
            // Find the bin edges. They are defined by the centers between the grid points + the outermost grid
            // points themselves.
            //       bin index, . is center, | is grid point
            //      000111111333
            // start|--.--|--.--|stop
            double leftEdge, rightEdge;
            if (bin == 0)
                leftEdge = frequencyv[0];
            else
                leftEdge = 0.5 * (frequencyv[bin - 1] + frequencyv[bin]);
            if (bin == frequencyv.size() - 1)
                rightEdge = frequencyv[bin];
            else
                rightEdge = 0.5 * (frequencyv[bin] + frequencyv[bin + 1]);

            // now integrate over all the points that fall between the bin edges, and divide by interval
            auto iLeft = TemplatedUtils::index(leftEdge, integrationGridv);
            auto iRight = TemplatedUtils::index(rightEdge, integrationGridv);
            double integral = TemplatedUtils::integrate<double>(integrationGridv, integrandv, iLeft, iRight);
            double average = factor * integral / (rightEdge - leftEdge);

            // Add this average to the correct point of the spectrum
            binnedSpectrumv[bin] += average;
        }
    }

    double LineProfile::integrateSpectrum(const Spectrum& spectrum, double spectrumMax) const
    {
        // Approximate by value at line center if line is much narrower than the resolution of
        // the provided SED, or if option forces this
        double spectrumResolutionAtCenter = spectrum.resolution(_center);
        if ((5 * _halfWidth_lorentz < spectrumResolutionAtCenter && 3 * _sigma_gauss < spectrumResolutionAtCenter)
            || Options::lineprofile_forceTrivialIntegration)
            return spectrum.evaluate(_center);

        // Else, do this very complicated integration with automated cutoff
        constexpr size_t numPoints = 20;
        const Array lineGrid = recommendedFrequencyGrid(numPoints);

        if (lineGrid.size() > 1 && Options::lineprofile_optimizedLineIntegration)
        {
            const size_t iCenter = numPoints / 2;
            // We also have the maximum of the spectrum to our disposal. As long as a frequency
            // interval * maxSpectrum is small, we can keep skipping points. Once a frequency
            // interval * maxSpectrum > than the threshold appears, re-evaluate the line profile at
            // this frequency point.

            // This value needs to be tuned to the desired accuracy. Some naive testing has shown
            // me that 1e-5 to 1e-6 is good for H2, while the H lines need a much smaller value of
            // 1e-9 to 1e-10.
            const double CUTOFFCONTRIBUTIONFRAC = 1e-9;
            double integral = 0;

            // We will use v1 and v2 outside of this function too
            double x2, v2, y2, x1, v1, y1;
            auto evaluateIntegralInterval = [&](size_t iRight) -> double {
                x2 = lineGrid[iRight];
                v2 = (*this)(x2);
                y2 = spectrum.evaluate(x2) * v2;

                x1 = lineGrid[iRight - 1];
                v1 = (*this)(x1);
                y1 = spectrum.evaluate(x1) * v1;

                return 0.5 * (x2 - x1) * (y2 + y1);
            };

            // Add the center contribution (iRight = iCenter, or [iCenter - 1, iCenter])
            if (iCenter > 0) integral += evaluateIntegralInterval(iCenter);

            // Integral over right wing (starting with [iCenter, iCenter + 1], or iRight = iCenter
            // + 1)
            for (size_t iRight = iCenter + 1; iRight < lineGrid.size(); iRight++)
            {
                double contribution = evaluateIntegralInterval(iRight);
                // cout << "interval " << iRight - 1 << ", " << iRight << '\n';
                integral += contribution;

                // Flag that breaks out of this loop if set to true by the end of this body
                bool stopIntegrating = false;

                double absContribution = abs(contribution);
                double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
                // If the contribution is insignificant
                if (absContribution < significanceThres)
                {
                    double last_voigt = max(v1, v2);

                    // Start skipping, until last_voigt * new_deltaX * sumSpec/2 is significant.
                    // Since we are moving further into the wing when skipping, we know that any
                    // new voig evaluation will be even smaller than the current one. So if v2 *
                    // deltaX * Ysum / 2 is insignificant, then the the same value but with the
                    // proper voigt function will be even smaller, and we can safely ignore
                    // everything until this value is big again (for example when a tall spectrum
                    // peak appears in the wing).
                    while (true)
                    {
                        // If the maximum * the interval all the way to the end * current voigt is
                        // insignificant, break off the calculation. (This value is much larger
                        // than the remaining contribution, so if even this value is insignificant
                        // we are safe to say that the rest of the integral does not matter)
                        double deltaNuToEnd = lineGrid[lineGrid.size() - 1] - lineGrid[iRight];

                        // Default value of spectrummax is 0 (this means don't use it)
                        if (spectrumMax && spectrumMax * last_voigt * deltaNuToEnd < significanceThres)
                        {
                            stopIntegrating = true;
                            break;
                        }

                        // The order of these statements is pretty important. I hope this works for
                        // all cases.
                        if (iRight >= lineGrid.size() - 1) break;
                        // The increment needs to happen right away, otherwise, when the loop is
                        // first entered, we make an estimate for a point that already has a value.
                        iRight++;
                        double deltaX = lineGrid[iRight] - lineGrid[iRight - 1];
                        double sumY = spectrum.evaluate(lineGrid[iRight]) + spectrum.evaluate(lineGrid[iRight - 1]);
                        double contributionGuess = last_voigt * deltaX * sumY * 0.5;
                        if (abs(contributionGuess) >= significanceThres)
                        {
                            iRight--;
                            break;
                        }
                    }
                }
                if (stopIntegrating) break;
            }

            // if iCenter is 1 or 0, then there is no left wing
            if (iCenter < 2) return integral;

            // Integral over left wing (loop stops after iRight == 1), starting with iRight =
            // iCenter - 1, or [iCenter - 2, iCenter - 1].
            for (size_t iRight = iCenter; iRight-- > 1;)
            {
                double contribution = evaluateIntegralInterval(iRight);
                // cout << "interval " << iRight - 1 << ", " << iRight << '\n';
                integral += contribution;

                bool stopIntegrating = false;

                double absContribution = abs(contribution);
                double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
                if (absContribution < significanceThres)
                {
                    double last_voigt = max(v1, v2);

                    double deltaNuToBegin = lineGrid[iRight] - lineGrid[0];
                    if (spectrumMax && last_voigt * spectrumMax * deltaNuToBegin < significanceThres)
                    {
                        stopIntegrating = true;
                        break;
                    }

                    while (true)
                    {
                        if (iRight <= 1) break;
                        iRight--;
                        double deltaX = lineGrid[iRight] - lineGrid[iRight - 1];
                        double sumY = spectrum.evaluate(lineGrid[iRight]) + spectrum.evaluate(lineGrid[iRight - 1]);
                        double contributionGuess = last_voigt * deltaX * sumY * 0.5;
                        if (abs(contributionGuess) >= significanceThres)
                        {
                            iRight++;
                            break;
                        }
                    }
                }
                if (stopIntegrating) break;
            }
            return integral;
        }  // Options::lineprofile_optimizedLineIntegration
        else
        {
            const Array& spectrumGrid = spectrum.frequencyv();
            Array frequencyv(spectrumGrid.size() + lineGrid.size());

            // Integration over the whole wavelength range. Merges the two grids, and writes the
            // result to the last argument
            merge(begin(spectrumGrid), end(spectrumGrid), begin(lineGrid), end(lineGrid), begin(frequencyv));

            Array integrandv(frequencyv.size());
            for (size_t i = 0; i < frequencyv.size(); i++)
            {
                double freq = frequencyv[i];
                double val = spectrum.evaluate(freq);
                // cout << val << '\n';
                integrandv[i] = val * (*this)(freq);
            }
            return TemplatedUtils::integrate<double>(frequencyv, integrandv);
        }
    }

    void LineProfile::simpleAddToBinned(const Array& frequencyv, Array& binnedSpectrumv, double factor) const
    {
        // find closest frequency
        auto right = upper_bound(begin(frequencyv), end(frequencyv), _center);

        if (right == begin(frequencyv) || right == end(frequencyv))
            // if out of range do nothing
            return;
        else
        {
            // else (index 1 to n-1)
            double dleft = _center - *(right - 1);
            double dright = *right - _center;

            // if right side is closest, this is the correct bin
            int i = distance(begin(frequencyv), right);
            // else, need the one to the left
            if (dleft < dright) i -= 1;

            binnedSpectrumv[i] += factor / (*right - *(right - 1));
        }
    }
}
