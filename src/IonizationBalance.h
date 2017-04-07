#ifndef _IONIZATIONBALANCE_H_
#define _IONIZATIONBALANCE_H_

#include "Array.h"
#include "Constants.h"

namespace Ionization
{
double ionizedFraction(double nH, double T, const Array& frequencyv,
                       const Array& specificIntensityv);

double crossSection(double frequency);

double recombinationRate(double T);

const double THRESHOLD = 3.28984196e15; // Hertz;
}

#endif /* _IONIZATIONBALANCE_H_ */
