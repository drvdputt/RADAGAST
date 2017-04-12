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

/* Heating due to thermalization of freed electrons */
double heating(double nH, double f, double T, const Array& frequency, const Array& specificIntensityv);

const double THRESHOLD = 3.28984196e15; // Hertz;
}

#endif /* _IONIZATIONBALANCE_H_ */
