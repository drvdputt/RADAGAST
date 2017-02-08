#include "NumUtils.h"
#include "IonizationBalance.h"

double Ionization::ionizedFraction(double nH, double T, vector<double> wavelength, vector<double> isrf)
{
    size_t Nlambda = wavelength.size();
    vector<double> integrand(Nlambda, 0);
    
    for (size_t n = 0; n < Nlambda && wavelength[n] < 912 * Constant::ANG_CM; n++)
    {
        integrand[n] = isrf[n] * wavelength[n] * crossSection(wavelength[n]);
    }
    double C =  recombinationRate(T) / Constant::PLANCK * NumUtils::integrate<double>(wavelength, integrand);
    
    // np^2 = (nH - np) * integral / alpha = (nH - np) * C
    // np^2 + C*np - C*nH = 0
    // np = (-C + sqrt(C^2 + 4 * C*nH)) / 2
    
    
    double np = (-C + sqrt(C*C + 4 * C * nH)) / 2.;
    return np / nH;
}

double Ionization::crossSection(double lambda)
{
    double E0 = 4.298e-1 / Constant::ERG_EV;
    double sigma0 = 5.475e4;
    double ya = 3.288e1;
    double P = 2.963;
    
    double x = Constant::PLANCKLIGHT / lambda / E0;
    double y = x;
    
    double Fy = (x - 1)*(x - 1) * pow(y, 0.5 * P - 5.5) * pow(1 + sqrt(y / ya), -P);
    
    return sigma0 * Fy * 1e-18;
}

double Ionization::recombinationRate(double T)
{
    double a = 7.982e-11;
    double b = 0.7480;
    double T0 = 3.148;
    double T1 = 7.036e5;
    
    return a / (sqrt(T/T0) * pow(1 + sqrt(T/T0), 1 - b) * pow(1 + sqrt(T/T1), 1 + b));
}

double Ionization::testIonization()
{
    // Wavelength grid to use for tests
    const size_t Nlambda = 10000;
    double lambdaMin = 0.0001 * Constant::UM_CM;
    double lambdaMax = 10 * Constant::UM_CM;
    
    vector<double> wavelength(Nlambda);
    double lambdaStepFactor = std::pow(lambdaMax / lambdaMin, 1. / Nlambda);
    float lambda = lambdaMin;
    for (size_t n = 0; n < Nlambda; n++)
    {
        wavelength[n] = lambda;
        lambda *= lambdaStepFactor;
    }
    
    // A blackbody
    double _bbTemp = 5.e4;
    vector<double> isrf = NumUtils::bbodyCGS<double>(wavelength, _bbTemp);
    
    // Convert to energy density and rescale
    for (double& d : isrf) d *= Constant::FPI / Constant::LIGHT * 1000000;
    
    return ionizedFraction(25, 500, wavelength, isrf);
}
