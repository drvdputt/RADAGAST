#ifndef _TESTING_H_
#define _TESTING_H_

#include "Array.h"

#include <vector>

namespace Testing
{
Array generateFrequencyGrid(size_t nFreq, double minFreq, double maxFreq);
Array freqToWavGrid(const Array& frequencyv);
void refineFrequencyGrid(Array& grid, size_t nPerLine, double spacingPower,
		std::vector<double> lineFreqv, std::vector<double> lineWidthv);
Array generateSpecificIntensity(const Array& frequencyv, double Tc, double G0);
Array freqToWavSpecificIntensity(const Array& frequencyv,
		const Array& specificIntensity_nu);

void testHydrogenCalculator();
void testIonizationCrossSection();

void testPhotoelectricHeating();
}

#endif /* _TESTING_H_ */
