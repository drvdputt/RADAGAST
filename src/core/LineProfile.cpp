#include "LineProfile.h"
#include "Constants.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"

#include <cmath>

using namespace std;

LineProfile::LineProfile(double center, double sigma_gauss, double halfWidth_lorenz)
                : _center{center}, _sigma_gauss{sigma_gauss},
                  _halfWidth_lorenz(halfWidth_lorenz)

{
	_one_sqrt2sigma = M_SQRT1_2 / _sigma_gauss;
	_a = _halfWidth_lorenz * _one_sqrt2sigma;
}

double LineProfile::operator()(double nu) const
{
	double x = (nu - _center) * _one_sqrt2sigma;
	// Note that the normalization factor is 1 / sqrt(2 pi sigma)
	return SpecialFunctions::voigt(_a, x) / Constant::SQRT2PI / _sigma_gauss;
}

Array LineProfile::recommendedFrequencyGrid(int numPoints, double width) const
{
	// Odd numbers of points make this easier
	if (!(numPoints % 2))
		numPoints++;

	// Because then this is the center
	int iCenter = numPoints / 2;

	// We will space the points according to a power law.
	double spacingPower = 2.5;
	// double total_distance = 1.e-6 * _center;
	double total_distance = width * (_halfWidth_lorenz + _sigma_gauss) / 2;
	double w = total_distance / pow(iCenter, spacingPower);

	// Fill in the values
	Array freqv(numPoints);
	freqv[iCenter] = _center;
	// Example: iCenter = 2
	// frequency index:
	// 0 1 2 3 4
	// i:
	// 2 1 0 1 2
	// So i must include iCenter --> use <=
	for (int i = 1; i <= iCenter; i++)
	{
		double d = w * pow(i, spacingPower);

		// This is safe because we forced the number of points to be odd
		freqv[iCenter + i] = _center + d;
		freqv[iCenter - i] = _center - d;
	}
	return freqv;
}

void LineProfile::addToBinned(const Array& frequencyv, Array& binnedSpectrumv,
                              double factor) const
{
	double gridMin = *begin(frequencyv);
	double gridMax = *(end(frequencyv) - 1);

	// Skip points of the line that fall outside of the grid
	const Array& lineGrid = recommendedFrequencyGrid();
	auto lineBegin = begin(lineGrid);
	auto lineEnd = end(lineGrid);
	while (*lineBegin < gridMin)
		lineBegin++;
	while (*(lineEnd - 1) > gridMax)
		lineEnd--;
	double lineMin = *lineBegin;
	double lineMax = *(lineEnd - 1);
	size_t numLinePoints = distance(lineBegin, lineEnd);

	// We will only do the calculation for frequency bins the line overlaps with

	// Find the left side of the bin for the first line point
	auto gridStart = upper_bound(begin(frequencyv), end(frequencyv), lineMin) - 1;
	// And the right side of the bin for the last line point (the above code guarantees that
	// lineMax is smaller than the max of frequencyv)
	auto gridStop = upper_bound(begin(frequencyv), end(frequencyv), lineMax);

	// The bin edges are defined by the centers between the grid points + the outermost grid
	// points themselves.
	//       bin index, . is center, | is grid point
	//      000111111333
	// start|--.--|--.--|stop
	size_t nCenters = distance(gridStart, gridStop);
	vector<double> binEdges;
	binEdges.reserve(nCenters);

	binEdges.emplace_back(*gridStart); // leftmost point
	for (auto right = gridStart + 1; right <= gridStop; right++)
	{
		auto left = right - 1;
		binEdges.emplace_back((*right + *left) / 2.); // centers
	}
	binEdges.emplace_back(*gridStop); // rightmost point

	// Now merge the bin edges and the line points, which will serve as an
	// integration grid
	Array integrationGridv(binEdges.size() + numLinePoints);
	merge(begin(binEdges), end(binEdges), lineBegin, lineEnd, begin(integrationGridv));

	// Calculate the integrand on this grid
	Array integrandv(integrationGridv.size());
	for (size_t i = 0; i < integrandv.size(); i++)
		integrandv[i] = (*this)(integrationGridv[i]);

	// Go over the bins, and integrate the parts of the line that fall within each bin

	// Start with the distance to the left edge of the first bin (the result will be added
	// to consecutive grid points 'offset', which, according to our algorithm, are located
	// at the centers of the bins.)
	size_t offset = distance(begin(frequencyv), gridStart);
	for (auto right = begin(binEdges) + 1; right != end(binEdges); right++, offset++)
	{
		// Find the correct integration range
		double leftBound = *(right - 1);
		double rightBound = *right;
		size_t iLeft = TemplatedUtils::index(leftBound, integrationGridv);
		size_t iRight = TemplatedUtils::index(rightBound, integrationGridv);

		// Integrate over this range
		double integral = TemplatedUtils::integrate<double>(integrationGridv,
		                                                    integrandv, iLeft, iRight);

		// Turn this integral into an average over the interval (= the bin)
		double average = factor * integral / (rightBound - leftBound);

		// Add this average to the correct point of the spectrum
		binnedSpectrumv[offset] += average;
	}
}

void LineProfile::addToSpectrum(const Array& frequencyv, Array& spectrumv, double factor) const
{
#define OPTIMIZED_LINE_ADD
#ifdef OPTIMIZED_LINE_ADD
	// Only calculate the voigt function for significant contributions of this line to the
	// total spectrum. Start at the line center (assumes there is a suitable frequency point
	// in the grid)
	size_t iCenter = TemplatedUtils::index(_center, frequencyv);
	const double CUTOFFWINGCONTRIBUTION = 1e-9;
	double wingThres = 1e-2 * factor * (*this)(frequencyv[iCenter]);

	// Add values for center and right wing:
	for (size_t i = iCenter; i < frequencyv.size(); i++)
	{
		double value = factor * (*this)(frequencyv[i]);
		spectrumv[i] += value;

		// Once we are in the wing (arbitrarily defined here), check if we can start
		// ignoring some points. Stop evaluating the wing if its contribution becomes
		// negligible
		double absval = abs(value);
		if (absval < wingThres && absval < abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
		{
			// Keep skipping points until the spectrum drops below the threshold for
			// significance of the line again, or until we reach the end of the
			// spectrum.
			while (true)
			{
				// If we reach te end of the spectrum, break. The main loop will
				// exit.
				if (i >= frequencyv.size() - 1)
					break;
				i++;
				// If the current value is significant compared to the spectrum
				// at index i, roll back by 1 and break. The main loop will then
				// increment back to the current i, and calculate the value of
				// the line profile for it.
				if (absval > abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
				{
					i--;
					break;
				}
				// If the current value is still too small compared to the
				// current spectrum at i (remember that our wing is descending
				// and becoming even smaller), skip this point.
			}
		}
	}

	if (iCenter == 0)
		return;

	// Add values for left wing. Stops when i == 0 (i-- decrements to i - 1, but returns the
	// original, so the loop safely stops when the body exits with i == 0)
	for (size_t i = iCenter; i-- > 0;)
	{
		double value = factor * (*this)(frequencyv[i]);
		spectrumv[i] += value;
		double absval = abs(value);
		// For an insignificant wing value
		if (value < wingThres && absval < abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
		{
			// Move through the spectrum until it the value of the spectrum is small
			// enough for the value of the wing point to be significant again.
			while (true)
			{
				if (i == 0)
					break;
				i--;
				if (absval > abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
				{
					// Found a significant value! Go back up, then
					// let the main loop continue.
					i++;
					break;
				}
				// Did not find significant value, move further in the
				// left wing.
			}
		}
	}
#else
	// Add the whole line to the spectrum
	for (size_t i = 0; i < frequencyv.size(); i++)
		spectrumv[i] += factor * (*this)(frequencyv[i]);
#endif /* OPTIMIZED_LINE_ADD */
}

double LineProfile::integrateSpectrum(const Spectrum& spectrum, double spectrumMax) const
{
	const Array& spectrumGrid = spectrum.frequencyv();
	const Array& lineGrid = recommendedFrequencyGrid(27, 5);
	Array frequencyv(spectrumGrid.size() + lineGrid.size());

	// Merges the two grids, and writes the result to the last argument
	merge(begin(spectrumGrid), end(spectrumGrid), begin(lineGrid), end(lineGrid),
	      begin(frequencyv));

	// Temporary override to check if the integration works with only the points generated
	// for the line
	frequencyv = lineGrid;

#define OPTIMIZED_LINE_INTEGRATION
#ifdef OPTIMIZED_LINE_INTEGRATION
	const size_t iCenter = TemplatedUtils::index(_center, frequencyv);
	// We also have the maximum of the spectrum to our disposal. As long as a frequency
	// interval * maxSpectrum is small, we can keep skipping points. Once a frequency
	// interval * maxSpectrum > than the threshold appears, re-evaluate the line profile at
	// this frequency point.

	// This value needs to be tuned to the desired accuracy. Some naive testing has shown me
	// that 1e-5 to 1e-6 is good for H2, while the H lines need a much smaller value of 1e-9
	// to 1e-10.
	const double CUTOFFCONTRIBUTIONFRAC = 1e-9;

	double integral = 0;

	double x2, v2, y2, x1, v1, y1;
	auto evaluateIntegralInterval = [&](size_t iRight) -> double {
		x2 = frequencyv[iRight];
		v2 = (*this)(x2);
		y2 = spectrum.evaluate(x2) * v2;

		x1 = frequencyv[iRight - 1];
		v1 = (*this)(x1);
		y1 = spectrum.evaluate(x1) * v1;

		return 0.5 * (x2 - x1) * (y2 + y1);
	};

	// Add the center contribution (iRight = iCenter, or [iCenter - 1, iCenter])
	if (iCenter > 0)
		integral += evaluateIntegralInterval(iCenter);

	// Integral over right wing (starting with [iCenter, iCenter + 1], or iRight = iCenter +
	// 1)
	for (size_t iRight = iCenter + 1; iRight < frequencyv.size(); iRight++)
	{
		double contribution = evaluateIntegralInterval(iRight);
		// cout << "interval " << iRight - 1 << ", " << iRight << endl;
		integral += contribution;

		// Flag that breaks out of this loop if set to true by the end of this body
		bool stopIntegrating = false;

		double absContribution = abs(contribution);
		double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
		// If the contribution is insignificant
		if (absContribution < significanceThres)
		{
			double last_voigt = max(v1, v2);

			/* Start skipping, until last_voigt * new_deltaX * sumSpec/2 is
			   significant. Since we are moving further into the wing when skipping,
			   we know that any new voig evaluation will be even smaller than the
			   current one. So if v2 * deltaX * Ysum / 2 is insignificant, then the
			   the same value but with the proper voigt function will be even
			   smaller, and we can safely ignore everything until this value is big
			   again (for example when a tall spectrum peak appears in the wing). */
			while (true)
			{
				/* If the maximum * the interval all the way to the end *
				   current voigt is insignificant, break off the calculation.
				   (This value is much larger than the remaining contribution,
				   so if even this value is insignificant we are safe to say
				   that the rest of the integral does not matter) */
				double deltaNuToEnd = frequencyv[frequencyv.size() - 1] -
				                      frequencyv[iRight];
				// Default value of spectrummax is 0 (this means don't use it)
				if (spectrumMax &&
				    spectrumMax * last_voigt * deltaNuToEnd < significanceThres)
				{
					stopIntegrating = true;
					break;
				}

				// The order of these statements is pretty important. I hope
				// this works for all cases.
				if (iRight >= frequencyv.size() - 1)
					break;
				// The increment needs to happen right away, otherwise, when the
				// loop is first entered, we make an estimate for a point that
				// already has a value.
				iRight++;
				double deltaX = frequencyv[iRight] - frequencyv[iRight - 1];
				double sumY = spectrum.evaluate(frequencyv[iRight]) +
				              spectrum.evaluate(frequencyv[iRight - 1]);
				double contributionGuess = last_voigt * deltaX * sumY * 0.5;
				if (abs(contributionGuess) >= significanceThres)
				{
					iRight--;
					break;
				}
			}
		}
		if (stopIntegrating)
			break;
	}

	// if iCenter is 1 or 0, then there is no left wing
	if (iCenter < 2)
		return integral;

	// Integral over left wing (loop stops after iRight == 1), starting with iRight =
	// iCenter - 1, or [iCenter - 2, iCenter - 1].
	for (size_t iRight = iCenter; iRight-- > 1;)
	{
		double contribution = evaluateIntegralInterval(iRight);
		// cout << "interval " << iRight - 1 << ", " << iRight << endl;
		integral += contribution;

		bool stopIntegrating = false;

		double absContribution = abs(contribution);
		double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
		if (absContribution < significanceThres)
		{
			double last_voigt = max(v1, v2);

			double deltaNuToBegin = frequencyv[iRight] - frequencyv[0];
			if (spectrumMax &&
			    last_voigt * spectrumMax * deltaNuToBegin < significanceThres)
			{
				stopIntegrating = true;
				break;
			}

			while (true)
			{
				if (iRight <= 1)
					break;
				iRight--;
				double deltaX = frequencyv[iRight] - frequencyv[iRight - 1];
				double sumY = spectrum.evaluate(frequencyv[iRight]) +
				              spectrum.evaluate(frequencyv[iRight - 1]);
				double contributionGuess = last_voigt * deltaX * sumY * 0.5;
				if (abs(contributionGuess) >= significanceThres)
				{
					iRight++;
					break;
				}
			}
		}
		if (stopIntegrating)
			break;
	}
	return integral;
#else
	Array integrandv(frequencyv.size());
	for (size_t i = 0; i < frequencyv.size(); i++)
	{
		double freq = frequencyv[i];
		double val = spectrum.evaluate(freq);
		// cout << val << endl;
		integrandv[i] = val * (*this)(freq);
	}
	return TemplatedUtils::integrate<double>(frequencyv, integrandv);
#endif /* OPTIMIZED_LINE_INTEGRATION */
}

namespace
{
LineProfile testLine()
{
	double c = 4.;
	double sg = 0.1;
	double hwl = 0.1;
	return LineProfile(c, sg, hwl);
}

Spectrum testSpectrum(double base)
{
	size_t numFreq = 100;
	Array frequencyv(numFreq);
	double fmin = .01;
	double fmax = 10.;
	double step = (fmax - fmin) / numFreq;
	for (size_t i = 0; i < numFreq; i++)
		frequencyv[i] = fmin + i * step;
	Array spectrumv(base, numFreq);
	return Spectrum(frequencyv, spectrumv);
}
} // namespace

void test_addToSpectrum()
{
	cout << "test_addToSpectrum" << endl;
	auto lp = testLine();

	double base = 1.;
	auto s = testSpectrum(base);
	Array frequencyv = s.frequencyv();
	Array spectrumv = s.valuev();

	double factor = 3.;
	lp.addToSpectrum(frequencyv, spectrumv, factor);

	double e = 1.e-15;
	for (size_t i = 0; i < frequencyv.size(); i++)
	{
		double linevalue = factor * lp(frequencyv[i]);
		double manual = base + linevalue;
		Error::fuzzyCheck("test value", spectrumv[i], manual, e);
	}
}

void test_integrateSpectrum()
{
	cout << "test_integrateSpectrum" << endl;
	auto lp = testLine();

	// Generate a flat spectrum of value base
	double base = 4.;
	auto s = testSpectrum(base);

	double integral = lp.integrateSpectrum(s);
	
	// This integral should be about equal to base, if gridpoints are chosen well
	double e = 1.e-3;
	Error::fuzzyCheck("test value", integral, base, e);
}
