#include "LineProfile.h"
#include "Constants.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"

#include <cmath>

using namespace std;

LineProfile::LineProfile(double center, double sigma_gauss, double halfWidth_lorenz)
                : _center{center}, _sigma_gauss{sigma_gauss}

{
	_one_sqrt2sigma = M_SQRT1_2 / _sigma_gauss;
	_a = halfWidth_lorenz * _one_sqrt2sigma;
}

double LineProfile::operator()(double nu) const
{
	double x = (nu - _center) * _one_sqrt2sigma;
	// Note that the normalization factor is 1 / sqrt(2 pi sigma)
	return SpecialFunctions::voigt(_a, x) / Constant::SQRT2PI / _sigma_gauss;
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
	double wingThres = 1e-6 * factor * (*this)(frequencyv[iCenter]);

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
				if (i >= frequencyv.size())
					break;
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
				i++;
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
			// Move through the spectrum until it the value of the spectrum
			// is small enough for the value of the wing point to be
			// significant again.
			while (true)
			{
				if (i == 0)
					break;
				if (absval > abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
				{
					// Found a significant value! Go back up, then
					// let the main loop continue.
					i++;
					break;
				}
				// Did not find significant value, move further in the
				// left wing.
				i--;
			}
		}
	}
#else
	// Add the whole line to the spectrum
	for (size_t i = 0; i < frequencyv.size(); i++)
		spectrumv[i] += factor * (*this)(frequencyv[i]);
#endif /* OPTIMIZED_LINE_ADD */
}

double LineProfile::integrateSpectrum(const Array& frequencyv, const Array& spectrumv) const
{
#define OPTIMIZED_LINE_INTEGRATION
#ifdef OPTIMIZED_LINE_INTEGRATION
	const size_t iCenter = TemplatedUtils::index(_center, frequencyv);
	// Get max of spectrum. As long as a frequency interval * maxSpectrum is small, we can
	// keep skipping points. Once a frequency interval * maxSpectrum > than the threshold
	// appears, re-evaluate the line profile at this frequency point.
	const double CUTOFFCONTRIBUTIONFRAC = 1e-9;

	double integral = 0;

	double x2, v2, y2, x1, v1, y1;
	auto evaluateIntegralInterval = [&](size_t iRight) -> double {
		x2 = frequencyv[iRight];
		v2 = (*this)(x2);
		y2 = spectrumv[iRight] * v2;

		x1 = frequencyv[iRight - 1];
		v1 = (*this)(x1);
		y1 = spectrumv[iRight - 1] * v1;

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
		double absContribution = abs(contribution);
		double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
		// If the contribution is insignificant
		if (absContribution < significanceThres)
		{
			double last_voigt = max(v1, v2);
			// The last voigt evaluation is v2 Start skipping, until last_voigt *
			// new_deltaX * sumSpec/2 is significant. Since we are moving further
			// into the wing when skipping, we know that any new voig evaluation
			// will be even smaller than the current one. So if v2 * deltaX * Ysum /
			// 2 is insignificant, then the the same value but with the proper voigt
			// function will be even smaller, and we can safely ignore everything
			// until this value is big again (for example when a tall spectrum peak
			// appears in the wing).
			while (true)
			{
				// The order of these statements is pretty important. I hope
				// this works for all cases.
				if (iRight >= frequencyv.size() - 1)
					break;
				// The increment needs to happen right away, otherwise, when the
				// loop is first entered, we make an estimate for a point that
				// already has a value.
				iRight++;
				double deltaX = frequencyv[iRight] - frequencyv[iRight - 1];
				double sumY = spectrumv[iRight] + spectrumv[iRight - 1];
				double contributionGuess = last_voigt * deltaX * sumY * 0.5;
				if (abs(contributionGuess) >= significanceThres)
				{
					iRight--;
					break;
				}
			}
		}
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
		double absContribution = abs(contribution);
		double significanceThres = abs(integral) * CUTOFFCONTRIBUTIONFRAC;
		if (absContribution < significanceThres)
		{
			double last_voigt = max(v1, v2);
			while (true)
			{
				if (iRight <= 1)
					break;
				iRight--;
				double deltaX = frequencyv[iRight] - frequencyv[iRight - 1];
				double sumY = spectrumv[iRight] - spectrumv[iRight - 1];
				double contributionGuess = last_voigt * deltaX * sumY * 0.5;
				if (abs(contributionGuess) >= significanceThres)
				{
					iRight++;
					break;
				}
			}
		}
	}

	return integral;
#else
	Array integrandv(spectrumv.size());
	for (size_t iRight = 0; iRight < frequencyv.size(); iRight++)
		integrandv[iRight] = spectrumv[iRight] * (*this)(frequencyv[iRight]);

	return TemplatedUtils::integrate<double>(frequencyv, integrandv);
#endif /* OPTIMIZED_LINE_INTEGRATION */
}
