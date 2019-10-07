#include "LineProfile.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "IOTools.hpp"
#include "SpecialFunctions.hpp"
#include "TemplatedUtils.hpp"

#include <cmath>
#include <iomanip>

using namespace std;

LineProfile::LineProfile(double center, double sigma_gauss, double halfWidth_lorentz)
                : _center{center}, _sigma_gauss{sigma_gauss},
                  _halfWidth_lorentz(halfWidth_lorentz)

{
	_one_sqrt2sigma = M_SQRT1_2 / _sigma_gauss;
	_a = _halfWidth_lorentz * _one_sqrt2sigma;
}

double LineProfile::operator()(double nu) const
{
	double x = (nu - _center) * _one_sqrt2sigma;
	// Note that the normalization factor is 1 / sqrt(2 pi) sigma
	return SpecialFunctions::voigt(_a, x) / Constant::SQRT2PI / _sigma_gauss;
}

Array LineProfile::recommendedFrequencyGrid(int numPoints) const
{
	// Odd numbers of points make this easier
	if (!(numPoints % 2))
		numPoints++;

	// Because then this is the center
	int iCenter = numPoints / 2;

	double fwhmGauss = 2.35 * _sigma_gauss;
	double fwhmLorentz = 2 * _halfWidth_lorentz;
	// Approximation from wikipedia (actual source was is 404)
	double fwhmVoigt = 0.5346 * fwhmLorentz + std::sqrt(0.2166 * fwhmLorentz * fwhmLorentz +
	                                                    fwhmGauss * fwhmGauss);

	double yMax = (*this)(_center);
	// With 2.33 sigma, we get about 0.99 of the gaussian
	double xMax = std::max(2.6 * _sigma_gauss,
	                       SpecialFunctions::lorentz_percentile(0.995, _halfWidth_lorentz));
	double yMin = (*this)(_center + xMax);
	double linearStep = (yMax - yMin) / (iCenter + 1); // +1 to avoid stepping below 0

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

		// Assume a gaussian. We want to go down a fixed step in the vertical direction.
		x = SpecialFunctions::inverse_gauss(y - linearStep, _sigma_gauss);

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

		if (x > fwhmVoigt / 2)
			break;
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
	// std::cout << "added " << left << " linear points" << std::endl;
	return freqv;
}

void LineProfile::addToBinned(const Array& frequencyv, Array& binnedSpectrumv,
                              double factor) const
{
	if (frequencyv.size() < 2)
		Error::runtime("frequencyv should have at least two values");

	double gridMin = frequencyv[0];
	double gridMax = frequencyv[frequencyv.size() - 1];

	const Array& lineGrid = recommendedFrequencyGrid();
	auto lineBegin = begin(lineGrid);
	auto lineEnd = end(lineGrid);

	// If line is completely outside grid, do nothing
	if (*lineBegin > gridMax || *(lineEnd - 1) < gridMin)
		return;

	// If some points for the line are outside the grid, ignore them
	while (*lineBegin < gridMin)
		lineBegin++;
	while (*(lineEnd - 1) > gridMax)
		lineEnd--;
	double lineMin = *lineBegin;
	double lineMax = *(lineEnd - 1);
	size_t numLinePoints = distance(lineBegin, lineEnd);

	// We will only do the calculation for frequency bins the line overlaps with.

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

double LineProfile::integrateSpectrum(const Spectrum& spectrum, double spectrumMax,
                                      std::string debug) const
{
	// Approximate by value at line center if line is much narrower than the resolution of
	// the provided SED
	double spectrumResolutionAtCenter = spectrum.resolution(_center);
	if (5 * _halfWidth_lorentz < spectrumResolutionAtCenter &&
	    3 * _sigma_gauss < spectrumResolutionAtCenter)
		return spectrum.evaluate(_center);

	// Else, do this very complicated integration with automated cutoff
	const Array& spectrumGrid = spectrum.frequencyv();
	const Array& lineGrid = recommendedFrequencyGrid(15);
	Array frequencyv(spectrumGrid.size() + lineGrid.size());

	// Merges the two grids, and writes the result to the last argument
	merge(begin(spectrumGrid), end(spectrumGrid), begin(lineGrid), end(lineGrid),
	      begin(frequencyv));

	// TODO: Temporary override
	frequencyv = lineGrid;

	if (!debug.empty())
	{
		auto f = IOTools::ofstreamFile(debug);
		for (double nu : frequencyv)
			f << std::setprecision(17) << nu << ' '
			  << SpecialFunctions::lorentz(nu - _center, _halfWidth_lorentz) << ' '
			  << SpecialFunctions::gauss(nu - _center, _sigma_gauss) << ' '
			  << (*this)(nu) << std::endl;
		f.close();
	}

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
