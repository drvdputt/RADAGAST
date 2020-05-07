// Make a wavelength grid based on the gas module settings (settings are fixed at compile time,
// for now)

#include "Constants.hpp"
#include "Functions.hpp"
#include "H2Data.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "SimpleColumnFile.hpp"

// parameters
const double lineWidthT = 20.;

const double min = 0.005 * Constant::UM;
const double max = 1000 * Constant::UM;

const double ppd_below70 = 100;
const double ppd_70to91 = 400;

// for H2 (91 to 120 nm)
const double ppd_91to120 = 200;

// for UV - IR (120 nm to 1 micron, for photoelectric heating amongst others)
const double ppd_120to1000 = 100;

// the rest of the IR / submm
const double ppd_IRtoSubmm = 50;

const std::vector<double> ppd_values = {ppd_below70, ppd_70to91, ppd_91to120, ppd_120to1000, ppd_IRtoSubmm};
const std::vector<double> ppd_right_bounds = {70 * Constant::NM, 91.2 * Constant::NM, 120 * Constant::NM,
                                              1 * Constant::UM, max};

// index of the ppd bound > lambda
int ppd_index(double lambda)
{
    auto iterator = std::upper_bound(std::begin(ppd_right_bounds), std::end(ppd_right_bounds), lambda);
    int index = std::distance(std::begin(ppd_right_bounds), iterator);
    return index;
}

// for line: at least a wavelength point on the line center, and 1 on each side to make sure the
// bin is small. There are already a couple of functions in the code that do something similar,
// but I want to redesign the algorithm a little.
void addLinePoints(double center, double width, std::vector<double>& output)
{
    // we want upper point to be center + 0.5 * width and lower point to be center - 0.5 * width
    // ~= center / (1 + 0.5 * width / center)
    double factor = (1 + 0.5 * width / center);
    double right = center * factor;
    double left = center / factor;
    output.push_back(left);
    output.push_back(center);
    output.push_back(right);
}

// Adds the point 'pointBelow' to the output vector, and a counterpoint (threshold * treshold /
// pointBelow) so that the geometric mean of pointBelow and counterpoint falls exactly on the
// threshold.
void addThresholdPoints(double threshold, double pointBelow, std::vector<double>& output)
{
    output.push_back(pointBelow);
    output.push_back(threshold * threshold / pointBelow);
}

void addContinuumPoints(double start, double stop, double ppd, std::vector<double>& output, bool skipStart,
                        bool skipStop)
{
    double ratio = stop / start;
    int numPoints = ceil(log10(ratio) * ppd);
    if (numPoints <= 0) return;
    double factor = pow(ratio, 1. / numPoints);

    if (!skipStart) output.push_back(start);

    double current = start;
    for (int i = 1; i < numPoints; i++)
    {
        current *= factor;
        output.push_back(current);
    }

    if (!skipStop) output.push_back(stop);
}

int main()
{
    // things to resolve (+ or - indicates priority): + H2 ion threshold (~ 80.3 nm, with
    // maximum at 72.3), ++ H ion threshold (91.2 nm), + H2 photodissociation continuum, ++ H2
    // photodissociation, - H2 rovibration lines, - H lines

    // algorithm: first, gather all the specialized frequency points and sort them. Then, fill
    // in any gaps using a certain amount of points per dex.
    std::vector<double> specialPoints;

    // H2 ion threshold (only used for opacity / radiation field)
    addThresholdPoints(80.3 * Constant::NM, 80.2 * Constant::NM, specialPoints);

    // H ion threshold (important for opacity / radiation field and H+/H balance)
    double threshold912 = Constant::LIGHT / GasModule::Ionization::THRESHOLD;
    addThresholdPoints(threshold912, 91.1 * Constant::NM, specialPoints);

    // H2 lines
    if (GasModule::Options::speciesmodelmanager_enableBigH2)
    {
        GasModule::H2Data h2Data(GasModule::Options::h2data_X_maxJ, GasModule::Options::h2data_X_maxV,
                                 GasModule::Options::h2data_E_maxJ, GasModule::Options::h2data_E_minV,
                                 GasModule::Options::h2data_E_maxV);
        int numLines;
        Array lineFreqv;
        Array lorentzWidthv;
        h2Data.lineInfo(numLines, lineFreqv, lorentzWidthv);
        int count = 0;
        for (int i = 0; i < numLines; i++)
        {
            double center = Constant::LIGHT / lineFreqv[i];
            if (center < 240 * Constant::NM)  // (UV only)
            {
                // for width, use formula copied from LevelCoefficients::lineProfile but with
                // constant temperature (200 K is a rough guildelinee).
                double sigma = center * GasModule::Functions::thermalVelocityWidth(lineWidthT, 2 * Constant::HMASS)
                               / Constant::LIGHT;
                addLinePoints(center, sigma, specialPoints);
                count++;
            }
        }
        std::cout << "Added points for " << count << " H2 (UV) lines\n";
    }

    std::sort(std::begin(specialPoints), std::end(specialPoints));

    std::vector<double> allpoints;

    // go over all the intervals between the specialized, and add log spaced points if necessary
    double start = min;
    for (int i = 0; i < specialPoints.size(); i++)
    {
        // determine current ppd index and quit if we're at the last threshold
        int ippd = ppd_index(start);
        if (ippd == ppd_right_bounds.size()) break;

        // get current ppd and next ppd boundary
        double ppd_at_start = ppd_values[ippd];
        double next_ppd_bound = ppd_right_bounds[ippd];

        // make sure we don't cross a ppd bound
        double stop = std::min(specialPoints[i], next_ppd_bound);

        // add the points with the ppd for this wavelength range
        allpoints.push_back(start);
        addContinuumPoints(start, stop, ppd_at_start, allpoints, true, true);

        // next iteration, start from the current endpoint
        start = stop;
    }
    // if we're not at the end yet, go over the remaining boundaries
    for (int ippd = ppd_index(start); ippd < ppd_right_bounds.size(); ippd++)
    {
        double stop = ppd_right_bounds[ippd];
        addContinuumPoints(start, stop, ppd_values[ippd], allpoints, false, true);
        start = stop;
    }
    // add the last point manually
    allpoints.push_back(max);

    // remove duplicates
    auto newEnd = std::unique(allpoints.begin(), allpoints.end());
    allpoints.resize(std::distance(allpoints.begin(), newEnd));

    // write out the result
    GasModule::OutColumnFile outfile("wlg.dat", {"wavelength (micron)", "frequency (hz)"}, 11);
    for (int i = 0; i < allpoints.size(); i++)
    {
        double lambda_micron = allpoints[i] / Constant::UM;
        double nu = Constant::LIGHT / allpoints[i];
        outfile.writeLine<std::vector<double>>({lambda_micron, nu});
    }
    return 0.;
}
