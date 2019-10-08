#ifndef CORE_OPTIONS_H_
#define CORE_OPTIONS_H_

namespace Options
{
// Write out the level matrices of h2 at setup (written to working_directory/h2/einsteinA.dat
// and levels.dat)
const bool h2fromfiles_plotLevelMatrices = false;

// When integrating the product of a line profile with a spectrum, use a cutoff heuristic to
// optimize this integration. When false, the line is integrated over the whole wavelength range
// instead.
const bool lineprofile_optimizedLineIntegration = true;

// When adding a line to a spectrum (i.e. averaging the line contribution over each bin of the
// spectrum), only evaluate the line profile where the contribution is significant. When false,
// the line is evaluated and added to the spectrum over the whole wavelength range.
const bool lineprofile_optimizedLineAdd = true;

} /* namespace Options */

#endif /* CORE_OPTIONS_H_ */
