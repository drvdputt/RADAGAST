#include "NumUtils.h"
#include "IonizationBalance.h"

double Ionization::ionizedFraction(double nH, double T, const vector<double>& wavelengthv, const vector<double>& isrf)
{
    size_t nWav = wavelengthv.size();
    vector<double> integrand(nWav, 0);

    for (size_t n = 0; n < nWav && wavelengthv[n] < 912 * Constant::ANG_CM; n++)
    {
        integrand[n] = isrf[n] * wavelengthv[n] * crossSection(wavelengthv[n]);
    }
    double C = NumUtils::integrate<double>(wavelengthv, integrand) / Constant::PLANCK / recombinationRate(T);

    // np^2 = (nH - np) * integral / alpha = (nH - np) * C
    // np^2 + C*np - C*nH = 0
    // np = (-C + sqrt(C^2 + 4 * C*nH)) / 2

    double np = (-C + sqrt(C*C + 4 * C * nH)) / 2.;
    return np / nH;
}

double Ionization::crossSection(double wavelength)
{
    double E0 = 4.298e-1 / Constant::ERG_EV;
    double sigma0 = 5.475e4;
    double ya = 3.288e1;
    double P = 2.963;

    double x = Constant::PLANCKLIGHT / wavelength / E0;
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
    const size_t nWav = 10000;
    double minWav = 0.0001 * Constant::UM_CM;
    double maxWav = 10 * Constant::UM_CM;

    vector<double> wavelengthv(nWav);
    double wavStepFactor = std::pow(maxWav / minWav, 1. / nWav);
    float wav = minWav;
    for (size_t n = 0; n < nWav; n++)
    {
        wavelengthv[n] = wav;
        wav *= wavStepFactor;
    }

    // A blackbody
    double _bbTemp = 5.e4;
    vector<double> isrf = NumUtils::bbodyCGS<double>(wavelengthv, _bbTemp);

    // Convert to energy density and rescale
    for (double& d : isrf) d *= Constant::FPI / Constant::LIGHT * 1000000;

    return ionizedFraction(25, 500, wavelengthv, isrf);
}
