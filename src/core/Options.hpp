#ifndef CORE_OPTIONS_H_
#define CORE_OPTIONS_H_

namespace Options
{
/////////////
// PHYSICS //
/////////////

// Take into account cooling of the gas by collisions with the grains. Also activates the
// corresponding heating term for the grains.
const bool cooling_gasGrainCollisions = true;

// To ensure self-consistency, add a source term to the upper level when the total recombination
// rate is larger than the sum of recombination rates to individual levels.
const bool hlevels_topoff = true;

// Turn on the 'recombination cooling' term from WD01
const bool grainphotoelectriceffect_recombinationCooling = false;

// Instead of using the Emin of WD01, use the expression provided in van Hoof (2004). Last time
// I tried this, weird things happened with the charge distribution, so please keep an eye on
// the charges if you decide to use this.
const bool weingartnerdraine2001_vanHoofEmin = false;

// Choose the number of electronic levels for H2. Give a number from 1 to 3 to include up to: 1
// B, 2 C+, 3 C-.
const int h2fromfiles_numExcitedLevels = 2;

//////////////////////////////////////////////
// NUMERIC / METHOD / ITERATION / PRECISION //
//////////////////////////////////////////////

// The maximum allowed amount of iterations for the (levels --> chemrates --> chemistry -->
// abundances --> levels)-loop.
const int densities_maxiterations = 25;

// When integrating the product of a line profile with a spectrum, use a cutoff heuristic to
// optimize this integration. When false, the line is integrated over the whole wavelength range
// instead.
const bool lineprofile_optimizedLineIntegration = true;

// When adding a line to a spectrum (i.e. averaging the line contribution over each bin of the
// spectrum), only evaluate the line profile where the contribution is significant. When false,
// the line is evaluated and added to the spectrum over the whole wavelength range.
const bool lineprofile_optimizedLineAdd = true;

///////////
// DEBUG //
///////////

// Print the system to solve when using the Eigen solver for the statistical equilibrium
// (currenly not used for H2).
const bool levelsolver_printEquations = false;

// Print the transition coefficients (A, B and C matrices) to standard output (only in DEBUG).
// Currenly quite cryptic, as it happens in the general class NLevel, and therefore the species
// is not mentioned in the output.
const bool nlevel_printLevelMatrices = false;

// VERY COSTLY: When integrating the B-coefficients, also do the integration 'manually' on a
// very fine grid, and report the difference between the optimized line integration and the
// 'manual' way. This is implemented in NLevel, but is actually a check on the quality of
// LineProfile::integrateSpectrum.
const bool nlevel_reportSpecIntegral = false;

// Write out the level matrices of h2 at setup (written to working_directory/h2/einsteinA.dat
// and levels.dat)
const bool h2fromfiles_plotLevelMatrices = false;

// Write out the loaded, singly interpolate, and doubly interpolated freeBound continuum data
// (gamma_nu) (written to working_directory/freebound/loadedContinuum.dat,
// T-interpolatedContinuum.dat, and gammanufb.dat).
const bool freebound_debugData = false;

// Write out the data loaded for the free-free Gaunt factor, both 2D and integrated over
// frequency (written to working_directory/freefree/gauntff.dat and integratedgauntff.dat)
const bool freefree_debugData = false;

} /* namespace Options */

#endif /* CORE_OPTIONS_H_ */
