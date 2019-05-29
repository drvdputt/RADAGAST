#ifndef _CONSTANTS_H
#define _CONSTANTS_H

namespace Constant
{
// "Physical" constants
//   updated for TRUST values 12 Jun 2015
constexpr double LIGHT = 2.99792458e10;       // Speed of light cm sec^-1
constexpr double PI = 3.14159265358979323846; // Pi
constexpr double BOLTZMAN = 1.3806488e-16;    // Boltzman Constant, erg K^-1
constexpr double PLANCK = 6.62606957e-27;     // Planck Constant erg sec
constexpr double SOL_RAD = 6.955e+10;         // Solar radius in cm.
constexpr double SOL_MASS = 1.9891e+33;       // Solar mass in gm.
constexpr double SOL_LUM = 3.839e33;           // Solar mass in erg / s
constexpr double HMASS_AMU = 1.0078250321;    // Hydrogen mass in atomic mass units.
constexpr double HMASS_CGS = 1.67353263e-24;  // Hydrogen mass in gm.

// Conversion factors and common constant combinations.
constexpr double HBAR = PLANCK / 2. / PI;
constexpr double PLANCKLIGHT = (PLANCK) * (LIGHT);   // hc in erg cm
constexpr double IPLANCKLIGHT = 1.0 / (PLANCKLIGHT); // 1/hc
constexpr double FPI = 4.0 * (PI);                   // 4pi
constexpr double FPI_HC = FPI * IPLANCKLIGHT;        // 4pi/hc
constexpr double FPI_HC2 = FPI_HC * IPLANCKLIGHT;    // 4pi/(hc)^2
constexpr double CGS_AMU = 6.02214151e+23;           // gm --> Atomic mass units
constexpr double AMU_CGS = 1.66053886e-24;           // Atomic mass units --> gm
constexpr double ANG_UM = 1.0e-4;                    // Angstrom --> um
constexpr double ANG_CM = 1.0e-8;                    // Angstrom --> cm
constexpr double UM_CM = 1.0e-4;                     // um --> cm
constexpr double CM_UM = 1.0e4;                      // cm --> um
constexpr double AU_CM = 1.49597871e+13;             // AU --> cm
constexpr double PC_CM = 3.08567758e+18;             // Parsec --> cm
constexpr double JY_CGS = 1.0e-23;                   // Constants without wave conversion.
constexpr double MKS_CGS = 10.0;
constexpr double U_J = (FPI) / (LIGHT);              // Convert Energy Density to Field
constexpr double VFAC = (FPI) / 3.0;                 // Volume factor (4/3)*pi
constexpr double ERG_EV = 6.24150974e11;             // Ergs -> electron volts

// New constants
// 1 Habing in ergs cm-3 (Energy density between 6 and 13.6 eV)
constexpr double HABING = 5.33e-14;
constexpr double ELECTRONMASS = 9.1094e-28;        // electron mass in grams
constexpr double PROTONMASS = 1.6726219e-24;       // proton mass in grams
constexpr double FINESTRUCTURE = .0072973525664;   // Fine structure constant
constexpr double RYDBERG = 13.605693009 / ERG_EV;  // Rydberg constant in ergs
// electron charge square in cgs (I think cgs uses Gaussian // units or something)
constexpr double ESQUARE = Constant::FINESTRUCTURE * Constant::PLANCKLIGHT / 2. /
			   Constant::PI;

constexpr double SQRT2PI = 2.50662827463; // sqrt(2*pi)
constexpr double CSQUARE_TWOPLANCK = Constant::LIGHT * Constant::LIGHT / 2. / Constant::PLANCK;
}

#endif
