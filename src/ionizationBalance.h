#ifndef _IONIZATIONBALANCE_H_
#define _IONIZATIONBALANCE_H_

#include<vector>

namespace Ionization 
{
    double ionizedFraction(double nH, double T, std::vector<double> wavelength, std::vector<double> isrf);
    
    double crossSection(double lambda);
    
    double recombinationRate(double T);
    
    double testIonization();
}

#endif
