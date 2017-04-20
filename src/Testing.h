#ifndef _TESTING_H_
#define _TESTING_H_

#include "Array.h"

#include <vector>

namespace Testing
{
std::vector<double> generateGeometricGridv(size_t nPoints, double min, double max);
std::vector<double> freqToWavGrid(const std::vector<double>& frequencyv);
void refineFrequencyGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower,
                         Array lineFreqv, Array lineWidthv);
Array generateSpecificIntensityv(const std::vector<double>& frequencyv, double Tc, double G0);
Array freqToWavSpecificIntensity(const std::vector<double>& frequencyv,
                                 const Array& specificIntensity_nu);

void testGasInterfaceImpl();
void testIonizationStuff();

void testPhotoelectricHeating();
}

#endif /* _TESTING_H_ */
