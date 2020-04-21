#include "FreeFree.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "TemplatedUtils.hpp"
#include <cassert>

using namespace std;

namespace
{
    // Some expected sizes of tables
    constexpr size_t NP_GAM2 = 81;
    constexpr size_t NP_U = 146;
    constexpr size_t NP_GAM2_INTEGRATED = 161;
    // See equations for free-free in gasphysics document, or in Rybicki and Lightman equations
    // 5.14a and 5.18a
    constexpr double e6 = Constant::ESQUARE * Constant::ESQUARE * Constant::ESQUARE;
    const double sqrt_2piOver3m = sqrt(2. * Constant::PI / 3. / Constant::ELECTRONMASS);
}  // namespace

namespace GasModule
{
    FreeFree::FreeFree()
    {
        readFullData();
        readIntegratedData();
    }

    /* Some of the code in this file (mainly the code for reading in the date) is an adaptation of
       the files interpolate3.c and interpolate4.c, which were downloaded from
       http://data.nublado.org/gauntff/ and included the Copyright statement below. */

    /*
      Copyright (c) 2014, Peter A.M. van Hoof.

      This program is provided 'as-is', without any express or implied warranty. In
      no event will the author be held liable for any damages arising from the use
      of this program.

      Permission is granted to anyone to use this program for any purpose, including
      commercial applications, and to alter it and redistribute it freely, subject
      to the following restrictions:

      1. The origin of this program must not be misrepresented; you must not claim
      that you created the original program. If you use this program in a product,
      an acknowledgment in the product documentation would be appreciated but
      is not required.
      2. Altered program versions must be plainly marked as such, and must not be
      misrepresented as being the original program.
      3. This notice may not be removed or altered from any further distribution.

      Peter A.M. van Hoof
      Royal Observatory of Belgium
      Ringlaan 3
      B-1180 Brussels
      Belgium
      p.vanhoof@oma.be
    */

    /*
      If you use any of these data in a scientific publication, please refer to

      van Hoof P.A.M., Williams R.J.R., Volk K., Chatzikos M., Ferland G.J., Lykins M., Porter R.L., Wang
      Y. Accurate determination of the free-free Gaunt factor, I -- non-relativistic Gaunt factors 2014,
      MNRAS, 444, 420

      van Hoof P.A.M., Ferland G.J., Williams R.J.R., Volk K., Chatzikos M., Lykins M., Porter R.L.
      Accurate determination of the free-free Gaunt factor, II -- relativistic Gaunt factors
      2015, MNRAS, 449, 2112
    */

    void FreeFree::readFullData()
    {
        /* Translated to c++ from interpolate3.c that came with the 2014 van Hoof paper (MNRAS 444
	   420) */
        ifstream input{IOTools::ifstreamRepoFile("dat/gauntff_merged_Z01.dat")};

        // buffer
        string line;

        /* skip copyright statement */
        while (getline(input, line))
        {
            if (line.at(0) != '#') break;
        }

        /* read magic number */
        const long gaunt_magic = 20140210L;
        const long gaunt_magic2 = 20140510L;
        long magic;
        istringstream(line) >> magic;
        if (magic != gaunt_magic && magic != gaunt_magic2)
            Error::runtime("read_table() found wrong magic number in file %s.\n");

        /* read dimensions of the table */
        size_t np_gam2, np_u;
        getline(input, line);
        istringstream(line) >> np_gam2 >> np_u;
        Error::equalCheck("Size and expected size of gam2 table", np_gam2, NP_GAM2);
        Error::equalCheck("Size and expected size of u table", np_u, NP_U);

        /* read start value for log(gamma^2) */
        getline(input, line);
        istringstream(line) >> _loggamma2Min;

        /* read start value for log(u) */
        getline(input, line);
        istringstream(line) >> _loguMin;

        /* read step size in dex */
        getline(input, line);
        istringstream(line) >> _logStep;

        _loggamma2Max = _loggamma2Min + static_cast<double>(np_gam2 - 1) * _logStep;
        _loguMax = _loguMin + static_cast<double>(np_u - 1) * _logStep;

        /* read atomic number when present */
        if (magic == gaunt_magic2) getline(input, line);

        /* next lines are comments */
        while (getline(input, line))
        {
            if (line.at(0) != '#') break;
        }

        /* the data */
        _fileGauntFactorvv.resize(np_u, np_gam2);
        for (size_t ipu = 0; ipu < np_u; ++ipu)
        {
            istringstream iss(line);
            for (size_t ipg2 = 0; ipg2 < np_gam2; ++ipg2)
            {
                double value;
                iss >> value;
                _fileGauntFactorvv(ipu, ipg2) = log(value);
            }
            getline(input, line);
        }
        DEBUG("Successfully read gauntff.dat" << endl);

        if (Options::freefree_debugData)
        {
            // Write out 2D guant factor data
            ofstream out = IOTools::ofstreamFile("freefree/gauntff.dat");
            for (size_t ipu = 0; ipu < np_u; ipu++)
            {
                for (size_t ipg2 = 0; ipg2 < np_gam2; ++ipg2) out << _fileGauntFactorvv(ipu, ipg2) << '\t';
                out << endl;
            }
            out.close();
        }
    }

    void FreeFree::readIntegratedData()
    {
        // Translated to c++ from interpolate3.c that came with the 2014 van Hoof paper (MNRAS 444
        // 420)
        ifstream input{IOTools::ifstreamRepoFile("dat/gauntff_freqint_Z01.dat")};

        // buffer
        string line;

        /* skip copyright statement */
        while (getline(input, line))
        {
            if (line.at(0) != '#') break;
        }

        /* read magic number */
        const long gaunt_magic = 20141008L;
        long magic;
        istringstream(line) >> magic;
        if (magic != gaunt_magic) Error::runtime("read_table() found wrong magic number in file %s.\n");

        /* read dimensions of the table */
        size_t np_gam2;
        getline(input, line);
        istringstream(line) >> np_gam2;
        Error::equalCheck("size and expected size of integrated gam2 table", np_gam2, NP_GAM2_INTEGRATED);

        /* read start value for log(gamma^2) */
        getline(input, line);
        istringstream(line) >> _loggamma2Min_integrated;

        /* read step size in dex */
        getline(input, line);
        istringstream(line) >> _logStep_integrated;
        _loggamma2Max_integrated = _loggamma2Min_integrated + static_cast<double>(np_gam2 - 1) * _logStep_integrated;

        /* read atomic number */
        getline(input, line);

        /* next lines are comments */
        while (getline(input, line))
        {
            if (line.at(0) != '#') break;
        }

        _fileGauntFactorv_integrated.resize(np_gam2);
        _loggamma2_integrated.resize(np_gam2);
        for (size_t ipg2 = 0; ipg2 < np_gam2; ++ipg2)
        {
            istringstream(line) >> _loggamma2_integrated[ipg2] >> _fileGauntFactorv_integrated[ipg2];
            _fileGauntFactorv_integrated[ipg2] = log(_fileGauntFactorv_integrated[ipg2]);
            getline(input, line);
        }

        if (Options::freefree_debugData)
        {
            ofstream out;
            out.open("freefree/integratedgauntff.dat");
            for (double logg2 = -5.9; logg2 < 9.9; logg2 += .1)
                out << logg2 << '\t' << integratedGauntFactor(logg2) << '\n';
            out.close();
        }
    }

    double FreeFree::gauntFactor(double logu, double logg2) const
    {
        // The data range should be more than enough. Return zero if out of range
        if (logu < _loguMin || _loguMax < logu || logg2 < _loggamma2Min || _loggamma2Max < logg2)
            return 0.;

        // Find the gamma^2-index to the right of logg2 , maximum the max column index)
        int iRight = ceil((logg2 - _loggamma2Min) / _logStep);
        // should be at least 1 (to extrapolate left)
        iRight = max(iRight, 1);
        // cannot be larger than the max column index
        iRight = min(iRight, static_cast<int>(_fileGauntFactorvv.size(1) - 1));
        int iLeft = iRight - 1;
        double xRight = _loggamma2Min + _logStep * iRight;
        double xLeft = xRight - _logStep;

        // Find the u-index above u
        int iUpper = ceil((logu - _loguMin) / _logStep);
        iUpper = max(iUpper, 1);
        iUpper = min(iUpper, static_cast<int>(_fileGauntFactorvv.size(0) - 1));
        int iLower = iUpper - 1;
        double yUp = _loguMin + _logStep * iUpper;
        double yLow = yUp - _logStep;

        double gff = TemplatedUtils::interpolateBilinear<double>(
            logg2, logu, xLeft, xRight, yLow, yUp, _fileGauntFactorvv(iLower, iLeft),
            _fileGauntFactorvv(iLower, iRight), _fileGauntFactorvv(iUpper, iLeft), _fileGauntFactorvv(iUpper, iRight));
        return exp(gff);
    }

    double FreeFree::integratedGauntFactor(double logg2) const
    {
        // Throw an error if out of range for now. Maybe allow extrapolation later.
        Error::rangeCheck("log(gamma^2)", logg2, _loggamma2Min_integrated, _loggamma2Max_integrated);

        /* Determine the indices of the data points we will interpolate between. Using the stored
           gamma^2 grid to interpolate using TemplatedUtils::evaluateLinInterpf would invoke a search
           algorithm for the nearest value. The fixed number of operations below should be faster. */

        // Find the gamma^2-index to the right of logg2 , maximum the max index)
        int iRight = ceil((logg2 - _loggamma2Min_integrated) / _logStep_integrated);
        // should be at least 1 (to extrapolate left)
        iRight = max(iRight, 1);
        // cannot be larger than the max column index
        iRight = min(iRight, static_cast<int>(_fileGauntFactorv_integrated.size() - 1));
        int iLeft = iRight - 1;

        double xRight = _loggamma2Min_integrated + _logStep_integrated * iRight;
        double xLeft = xRight - _logStep_integrated;
        double fLeft = _fileGauntFactorv_integrated[iLeft];
        double fRight = _fileGauntFactorv_integrated[iRight];
        double loggff_integrated_ip = TemplatedUtils::interpolateLinear(logg2, xLeft, xRight, fLeft, fRight);
        return exp(loggff_integrated_ip);
    }

    void FreeFree::addEmissionCoefficientv(double T, const Array& eFrequencyv, Array& gamma_nuv) const
    {
        // 32pi e^6 / 3mc^3 * sqrt(2pi / 3m)
        constexpr double c3 = Constant::LIGHT * Constant::LIGHT * Constant::LIGHT;
        const double gamma_nu_constantFactor =
            32 * Constant::PI * e6 / 3. / Constant::ELECTRONMASS / c3 * sqrt_2piOver3m;

        // gamma is fixed for a given temperature
        double kT = Constant::BOLTZMAN * T;
        double sqrtkT = sqrt(kT);
        double logg2 = log10(Constant::RYDBERG / kT);

        // this emissivity is smooth enough to just pick the representative frequency in each bin
        for (size_t iFreq = 0; iFreq < eFrequencyv.size(); iFreq++)
        {
            double u = Constant::PLANCK * eFrequencyv[iFreq] / kT;
            double logu = log10(u);
            double gammaNu = gamma_nu_constantFactor / sqrtkT * exp(-u) * gauntFactor(logu, logg2);
            gamma_nuv[iFreq] += gammaNu;
        }
    }

    void FreeFree::addOpacityCoefficientv(double T, const Array& oFrequencyv, Array& opCoeffv) const
    {
        for (size_t iFreq = 0; iFreq < oFrequencyv.size(); iFreq++)
        {
            double nu = oFrequencyv[iFreq];
            opCoeffv[iFreq] += opacityCoefficient(nu, T);
        }
    }

    double FreeFree::opacityCoefficient(double nu, double T) const
    {
        // C = 4 e^6 / 3mhc * sqrt(2pi / 3m)
        const double opCoef_constantFactor =
            4 * e6 / 3. / Constant::ELECTRONMASS / Constant::PLANCK / Constant::LIGHT * sqrt_2piOver3m;
        double kT = Constant::BOLTZMAN * T;
        double sqrtkT = sqrt(kT);
        double loggamma2 = log10(Constant::RYDBERG / kT);
        double u = Constant::PLANCK * nu / kT;
        double logu = log10(u);

        // C / nu^3 (1 - exp(-u)) gff(u, gamma^2)
        double opCoeffNu = opCoef_constantFactor / sqrtkT / nu / nu / nu * -expm1(-u) * gauntFactor(logu, loggamma2);
        return opCoeffNu;
    }

    double FreeFree::cooling(double np_ne, double T) const
    {
        // safety for log
        if (T <= 0) return 0;

        // Constant factor from 1998-Sutherland equation 18
        constexpr double fk = 1.42554e-27;
        // Interpolate and use the frequency integrated gaunt factor
        double logg2 = log10(Constant::RYDBERG / Constant::BOLTZMAN / T);
        double from_data = np_ne * fk * sqrt(T) * integratedGauntFactor(logg2);
        return from_data;
    }
}
