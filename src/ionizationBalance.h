#ifndef _IONIZATIONBALANCE_H_
#define _IONIZATIONBALANCE_H_

#include "Constants.h"

#include <vector>

namespace Ionization
{
double ionizedFraction(double nH, double T, const std::vector<double>& frequencyv,
		const std::vector<double>& specificIntensityv);

double crossSection(double frequency);

double recombinationRate(double T);

const double ionizationThreshold = Constant::LIGHT / 912 / Constant::ANG_CM;
}

#endif
