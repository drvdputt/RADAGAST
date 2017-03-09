#ifndef _TESTING_H_
#define _TESTING_H_
#include <vector>

namespace Testing
{
std::vector<double> generateFrequencyGrid(size_t nFreq, double minFreq, double maxFreq);
std::vector<double> freqToWavGrid(const std::vector<double>& frequencyv);
void refineFrequencyGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower,
		std::vector<double> lineFreqv, std::vector<double> lineWidthv);
std::vector<double> generateSpecificIntensity(const std::vector<double>& frequencyv, double Tc, double G0);
std::vector<double> freqToWavSpecificIntensity(const std::vector<double>& frequencyv,
		const std::vector<double>& specificIntensity_nu);

void testGasSpecies();
void testIonizationCrossSection();

void testPhotoelectricHeating();
}

#endif
