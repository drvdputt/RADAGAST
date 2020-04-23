#ifndef CORE_TESTING_HPP
#define CORE_TESTING_HPP

#include "Array.hpp"
#include "Constants.hpp"
#include "GasInterface.hpp"
#include <string>
#include <vector>

namespace GasModule
{
    class FreeBound;
    class GasInterfaceImpl;
    class GrainSolution;
    class LevelCoefficients;
    class Spectrum;

    namespace Testing
    {
        // ABSORPTION EFFICIENCIES FROM FILE //
        std::vector<Array> qAbsvvForTesting(bool car, const Array& av, const Array& frequencyv);

        // Makeshift dust equilibrium calculation
        double equilibriumTemperature(const Array& frequencyv, const Array& meanIntensityv, const Array& crossSectionv,
                                      double minT, double maxT);

        // USEFUL CONSTANTS //
        constexpr double defaultMinFreq = Constant::LIGHT / (1e3 * Constant::UM);
        constexpr double defaultMaxFreq = Constant::LIGHT / (0.005 * Constant::UM);

        // UTILITY FUNCTIONS //
        Array generateGeometricGridv(size_t nPoints = 200, double min = 1e11, double max = 1e16);

        std::vector<double> freqToWavGrid(const std::vector<double>& frequencyv);

        Array defaultFrequencyv(size_t numPoints = 100);

        /** Inserts frequency points into a vector of frequencies, given a number of points per
            line, a power for the subgrid per line, and a center and characteristic with per
            line. */
        void refineFrequencyGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower, Array lineFreqv,
                                 Array freqWidthv);

        /** This clumsy thing should give us a grid with some extra points in the correct
            locations. */
        Array improveFrequencyGrid(const LevelCoefficients& boundBound, const Array& oldpoints);
        Array improveFrequencyGrid(const FreeBound& freeBound, const Array& oldPoints);

        // TEST RUNS //
        /** Calculates an equilibrium GasState using the given GasInterface and writes out some
            results in the given directory (relative to the working directory). The environment can
            be specified through the optional arguments. The ambient radiation field has a color
            temperature of Tc and a mean UV intensity of G0 habing, and the gas density is n. The
            initial temperature can also be chosen. */
        void runGasInterfaceImpl(const GasModule::GasInterface& gi, const std::string& outputPath, double Tc = 20000,
                                 double G0 = 2, double n = 10);

        void writeGasState(const std::string& outputPath, const GasModule::GasInterface& gi,
                           const GasModule::GasState& gs);

        /** Write out the grain properties to a file. One file will be written per grain
            population, containing one line per size bin. The columns are the size, number density,
            mass density and temperature. The last argument @c bulkCar indicates which bulk mass
            density to use for converting the number density for each bin into a mass density. */
        void writeGrains(const std::string& outputPath, const std::vector<GrainSolution>& grs, bool bulkCar);

        void plotHeatingCurve_main();
        void plotHeatingCurve(const GasModule::GasInterface& gi, const std::string& outputPath, double n,
                              const Spectrum& meanIntensity, GasModule::GrainInterface&);

        void plotIonizationStuff();

        /** Write out data files to recreate figure 2 of PS64. Current results: the curve for n=5
            looks good, but the last point (l=3) of the n=4 curve was a bit to high. When
            retrieving the PS64 coefficients, a recursive relation is used, and one can choose to
            start from either l = 0 or l = n - 1. By applying both methods, and taking the average
            of the two results, the figure is approached quite well. */
        void plotPS64Collisions();

        /** Try to recreate the efficiency plot of WD01. Writes output to $(pwd)/photoElectric/ */
        void plotPhotoelectricHeating();

        /** Interpolate some simple functions, and plot the results. Used for eyeballing
            correctness of several interpolation routines. */
        void plotInterpolationTests();

        void plotChemistryTest();

        /** Test the H2 implementation. */
        void runH2(bool write);

        /** Generates a GasInterface object with the maximum number of levels and a suitable
            frequency grid. */
        GasModule::GasInterface genFullModel(bool refine = false);

        GasModule::GasInterface genHonlyModel();

        void genMRNDust(GasModule::GrainInterface&, double nHtotal, const Spectrum& meanIntensityv, bool car);

        /** Run the model with the maximum number of levels, and iterate until the temperature
            converges. */
        void runFullModel();

        /** A test with many dust grain sizes, following the MRN distribution (see Weingartner
            2001, ApJ, 548, 296). */
        void runMRNDust(bool write = true, double nH = 1.e4, double Tc = 1.e4, double lumSol = 1.e2,
                        bool own_dir = false);
    }  // namespace Testing
}

#endif  // CORE_TESTING_HPP
