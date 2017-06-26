#ifndef _TESTING_H_
#define _TESTING_H_

#include "Array.h"

#include <vector>

namespace Testing
{
// UTILITY FUNCTIONS //
std::vector<double> generateGeometricGridv(size_t nPoints, double min, double max);
std::vector<double> freqToWavGrid(const std::vector<double>& frequencyv);
void refineFrequencyGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower,
                         Array lineFreqv, Array lineWidthv);
Array generateSpecificIntensityv(const std::vector<double>& frequencyv, double Tc, double G0);
Array freqToWavSpecificIntensity(const std::vector<double>& frequencyv,
                                 const Array& specificIntensity_nu);

// TEST RUNS //
void testGasInterfaceImpl();
void testIonizationStuff();

/* Try to recreate the efficiency plot of WD01. Writes output to $(pwd)/photoElectric/*/
void testPhotoelectricHeating();

/* Writes out the A-coefficients of a fully collapsed H-model, so that they can be compared to the
   NIST values between different n. */
void testACollapse();

/* Write out data files to recreate figure 2 of PS64. Current results: the curve for n=5 looks good,
   but the last point (l=3) of the n=4 curve is a bit to high. This might be because the A
   coefficients used by Pengelley were different. I can't track them down easily though, as he
   refers to another paper, which in turn refers to his thesis. */
void testPS64Collisions();

/* Perform a direct comparison of the two LevelDataProvider classes HydrogenHardcoded and
   HydrogenFromFiles. */
void compareFromFilesvsHardCoded();
}

#endif /* _TESTING_H_ */
