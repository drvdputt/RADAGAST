#ifndef _IONIZATIONBALANCE_H_
#define _IONIZATIONBALANCE_H_

#include<vector>

namespace Ionization
{
    double ionizedFraction(double nH, double T, const std::vector<double>& wavelength, const std::vector<double>& isrf);

    double crossSection(double lambda);

    double recombinationRate(double T);

    double testIonization();
}

#endif
