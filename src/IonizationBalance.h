#ifndef _IONIZATIONBALANCE_H_
#define _IONIZATIONBALANCE_H_

#include "Array.h"

namespace Ionization
{
/** Solves the ionization balance (in the nebular approximation, i.e. this function assumes that all
    hydrogen is in the ground state). */
double solveBalance(double nH, double T, const Array& frequencyv, const Array& specificIntensityv);

/** Photoionization rate coefficient. Integrates over the spectrum. Multiply with neutral density
    to get total rate. */
double photoRateCoeff(const Array& frequencyv, const Array& specificIntensityv);

/** Photoionization cross section in cm2 */
double crossSection(double frequency);

/** Collisional ionization rate coefficient from the ground state in cm3 / s */
double collisionalRateCoeff(double T);

/** Total radiative recombination rate coefficient in cm3 / s */
double recombinationRateCoeff(double T);

/** Heating due to thermalization of freed electrons */
double heating(double nH, double f, double T, const Array& frequencyv,
               const Array& specificIntensityv);

/** Kinetic energy lost during recombination (basically the photon energy, minus the binding energy
    contribution. */
double cooling(double nH, double f, double T);

const double THRESHOLD = 3.28984196e15; // Hertz;
}

#endif /* _IONIZATIONBALANCE_H_ */
