#ifndef _CONSTANTS_H
#define _CONSTANTS_H

namespace Constant
{
    // "Physical" constants
    //   updated for TRUST values 12 Jun 2015
    constexpr double LIGHT = 2.99792458e10;        // Speed of light cm sec^-1
    constexpr double PI = 3.14159265358979323846;  // Pi
    constexpr double BOLTZMAN = 1.3806488e-16;     // Boltzman Constant, erg K^-1
    constexpr double PLANCK = 6.62606957e-27;      // Planck Constant erg sec
    constexpr double SOL_RAD = 6.955e+10;          // Solar radius in cm.
    constexpr double SOL_MASS = 1.9891e+33;        // Solar mass in gm.
    constexpr double SOL_LUM = 3.839e33;           // Solar mass in erg / s
    constexpr double HMASS_AMU = 1.0078250321;     // Hydrogen mass in atomic mass units.
    constexpr double HMASS = 1.67353263e-24;       // Hydrogen mass in gm.

    // Conversion factors and common constant combinations.
    constexpr double HBAR = PLANCK / 2. / PI;
    constexpr double PLANCKLIGHT = PLANCK * LIGHT;  // hc in erg cm
    constexpr double FPI = 4.0 * PI;                // 4pi
    constexpr double AMU = 1.66053886e-24;          // Atomic mass units --> gm
    constexpr double ANGSTROM = 1.0e-8;             // 1 Angstrom in cm
    constexpr double UM = 1.0e-4;                   // 1 micrometer in cm
    constexpr double AU = 1.49597871e+13;           // AU --> cm
    constexpr double PC = 3.08567758e+18;           // Parsec --> cm
    constexpr double ERG_EV = 6.24150974e11;        // Ergs -> electron volts
    constexpr double EV = 1.60218e-12;             // electron volt in erg

    // New constants
    // 1 Habing in ergs cm-3 (Energy density between 6 and 13.6 eV)
    constexpr double HABING_FLUX = 1.6e-3;             // erg s-1 cm-2
    constexpr double HABING_DENS = 5.33e-14;           // erg s-1 cm-3, divide the above by c
    constexpr double ELECTRONMASS = 9.1094e-28;        // electron mass in grams
    constexpr double PROTONMASS = 1.6726219e-24;       // proton mass in grams
    constexpr double FINESTRUCTURE = .0072973525664;   // Fine structure constant
    constexpr double RYDBERG = 13.605693009 / ERG_EV;  // Rydberg constant in ergs
    // electron charge square in cgs (I think cgs uses Gaussian // units or something)
    constexpr double ESQUARE = Constant::FINESTRUCTURE * Constant::PLANCKLIGHT / 2. / Constant::PI;

    constexpr double SQRT2PI = 2.50662827463;  // sqrt(2*pi)
    constexpr double CSQUARE_TWOPLANCK = Constant::LIGHT * Constant::LIGHT / 2. / Constant::PLANCK;
}

#endif
