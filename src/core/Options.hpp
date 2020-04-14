#ifndef CORE_OPTIONS_HPP
#define CORE_OPTIONS_HPP

namespace GasModule
{
    namespace Options
    {
        /////////////
        // PHYSICS //
        /////////////

        const bool speciesmodelmanager_enableBigH2 = true;

        // Determine which levels and transition data to load/use for the electronic ground
        // state (X) and the electronically excited states (E). These can be set separately
        // because different levels are important for different processes.
        const int h2data_X_maxJ = 6;
        const int h2data_X_maxV = 1;
        const int h2data_E_maxJ = 4;
        const int h2data_E_minV = 8;
        const int h2data_E_maxV = 25;

        // Take into account cooling of the gas by collisions with the grains. Also activates the
        // corresponding heating term for the grains.
        const bool cooling_gasGrainCollisions = true;

        // To ensure self-consistency, add a source term to the upper level when the total
        // recombination rate is larger than the sum of recombination rates to individual levels.
        const bool hlevels_topoff = true;

        // Turn on the 'recombination cooling' term from WD01. Off by default, because I believe
        // this might count some things double in the equations for gas-grain collisional cooling.
        const bool grainphotoelectriceffect_recombinationCooling = false;

        // Instead of using the Emin of WD01, use the expression provided in van Hoof (2004). The
        // first time I tried this, weird things happened with the charge distribution. Things look
        // fine now, but keep an eye on the average charge.
        const bool weingartnerdraine2001_vanHoofEmin = true;

        // Choose the number of electronic levels for H2. Give a number from 1 to 3 to include up
        // to: 1 B, 2 C+, 3 C-.
        const int h2data_numExcitedLevels = 2;

        const int grainsolution_minGrainTemp = 1.;
        const int grainsolution_maxGrainTemp = 300.;
        const int grainsolution_tPointsPerDex = 100.;
        const int grainphotoelectriccalculator_maxcharges = 16;

        //////////////////////////////////////////////
        // NUMERIC / METHOD / ITERATION / PRECISION //
        //////////////////////////////////////////////

        // The maximum allowed amount of iterations for the (levels --> chemrates --> chemistry -->
        // abundances --> levels)-loop.
        const int densities_maxiterations = 25;

        // When integrating the product of a line profile with a spectrum, use a cutoff heuristic
        // to optimize this integration. When false, the line is integrated over the whole
        // wavelength range instead.
        const bool lineprofile_optimizedLineIntegration = true;

        // When adding a line to a spectrum (i.e. averaging the line contribution over each bin of
        // the spectrum), only evaluate the line profile where the contribution is significant.
        // When false, the line is evaluated and added to the spectrum over the whole wavelength
        // range.
        const bool lineprofile_optimizedLineAdd = true;

        ///////////
        // DEBUG //
        ///////////

        // Print the system to solve when using the Eigen solver for the statistical equilibrium
        // (currenly not used for H2).
        const bool levelsolver_printEquations = false;

        // Print the transition coefficients (A, B and C matrices) to standard output (only in
        // DEBUG). Currenly quite cryptic, as it happens in the general class, and therefore the
        // species is not mentioned in the output.
        const bool levelcoefficients_printLevelMatrices = false;

        // VERY COSTLY: When integrating the B-coefficients, also do the integration 'manually' on
        // a very fine grid, and report the difference between the optimized line integration and
        // the 'manual' way. This is implemented in LevelCoefficients, but is actually a check on
        // the quality of LineProfile::integrateSpectrum.
        const bool levelcoefficients_reportSpecIntegral = false;

        // Write out the level matrices of h2 at setup (written to
        // working_directory/h2/einsteinA.dat and levels.dat)
        const bool h2data_plotLevelMatrices = false;

        // Write out the loaded, singly interpolate, and doubly interpolated freeBound continuum
        // data (gamma_nu) (written to working_directory/freebound/loadedContinuum.dat,
        // T-interpolatedContinuum.dat, and gammanufb.dat).
        const bool freebound_debugData = false;

        // Write out the data loaded for the free-free Gaunt factor, both 2D and integrated over
        // frequency (written to working_directory/freefree/gauntff.dat and integratedgauntff.dat)
        const bool freefree_debugData = false;

    } /* namespace Options */
}
#endif  // CORE_OPTIONS_HPP
