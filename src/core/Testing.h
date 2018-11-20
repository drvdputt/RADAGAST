#ifndef GASMODULE_GIT_SRC_TESTING_H_
#define GASMODULE_GIT_SRC_TESTING_H_

#include "Array.h"
#include "Constants.h"
#include "GasInterface.h"

#include <string>
#include <vector>

class FreeBound;
class GasInterfaceImpl;
class NLevel;
class Spectrum;

namespace Testing
{
// ABSORPTION EFFICIENCY FROM FILE //
/* Some static hacks. The first one reads some data ripped from SKIRT into variables declared in
   the cpp file. Then, a Qabsv for a single size can be generated using generateQabsv. The third
   function, qAbsvvForTesting does the read step, and then calls the second function for every
   size of the list. */
void readQabs(bool car);

Array generateQabsv(double a, const Array& frequencyv);

std::vector<Array> qAbsvvForTesting(const Array& av, const Array& frequencyv);

// USEFUL CONSTANTS //
constexpr double defaultMinFreq = Constant::LIGHT / (1e3 * Constant::UM_CM);
constexpr double defaultMaxFreq = Constant::LIGHT / (0.005 * Constant::UM_CM);

// UTILITY FUNCTIONS //

Array generateGeometricGridv(size_t nPoints = 200, double min = 1e11, double max = 1e16);

std::vector<double> freqToWavGrid(const std::vector<double>& frequencyv);

Array freqToWavGrid(const Array& frequencyv);

Array defaultCoarseFrequencyv();

/** Inserts frequency points into a vector of frequencies, given a number of points per line, a
    power for the subgrid per line, and a center and characteristic with per line. */
void refineFrequencyGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower,
                         Array lineFreqv, Array freqWidthv);

Array generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0);

Array freqToWavSpecificIntensity(const Array& frequencyv, const Array& specificIntensity_nu);

/** This clumsy thing should give us a grid with some extra points in the correct locations. */
Array improveFrequencyGrid(const NLevel& boundBound, const Array& oldpoints);
Array improveFrequencyGrid(const FreeBound& freeBound, const Array& oldPoints);

// TEST RUNS //
/** Calculates an equilibrium GasState using the given GasInterface and writes out some results
    in the given directory (relative to the working directory). The environment can be specified
    through the optional arguments. The ambient radiation field has a color temperature of Tc
    and a mean UV intensity of G0 habing, and the gas density is n. The initial temperature can
    also be chosen. */
void runGasInterfaceImpl(const GasModule::GasInterface& gi, const std::string& outputPath,
                         double Tc = 20000, double G0 = 2, double n = 10,
                         double expectedTemperature = 8000);

void writeGasState(const std::string& outputPath, const GasModule::GasInterface& gi,
                   const GasModule::GasState& gs);

void plotHeatingCurve_main();
void plotHeatingCurve(const GasInterfaceImpl& gi, const std::string& outputPath,
                      const Spectrum& specificIntensity, double n);

void plotIonizationStuff();

/** Write out data files to recreate figure 2 of PS64. Current results: the curve for n=5 looks
    good, but the last point (l=3) of the n=4 curve was a bit to high. When retrieving the PS64
    coefficients, a recursive relation is used, and one can choose to start from either l = 0 or
    l = n - 1. By applying both methods, and taking the average of the two results, the figure
    is approached quite well. */
void plotPS64Collisions();

/** Try to recreate the efficiency plot of WD01. Writes output to $(pwd)/photoElectric/ */
void plotPhotoelectricHeating();

/** Interpolate some simple functions, and plot the results. Used for eyeballing correctness of
    several interpolation routines. */
void plotInterpolationTests();

/** Test the H2 implementation. */
void runH2(bool write);

/** Do two full runs (determine 1 equilibrium GasState) using both HFF and HHC. The results will
    be written out to separate files. */
void runFromFilesvsHardCoded();

/** Generates a GasInterface object with the maximum number of levels and a suitable frequency
    grid. */
GasModule::GasInterface genFullModel();

GasModule::GasInterface genHonlyModel();

/** Run the model with the maximum number of levels, and iterate until the temperature
    converges. */
void runFullModel();

/** A test run for the model with a fixed dust model included. */
void runWithDust(bool write = true);

void plotHlines();

void plotH2lines();
} // namespace Testing

#endif /* GASMODULE_GIT_SRC_TESTING_H_ */