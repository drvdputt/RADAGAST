/* Some of the code in this file is adapted from the files interpolate3.c and interpolate4.c, which
 were downloaded from http://data.nublado.org/gauntff/ and included the Copyright statement below.
 */

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

#include "FreeFree.h"
#include "Constants.h"
#include "Error.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "Testing.h"
#include "global.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

#define NP_GAM2 81
#define NP_U 146

#define NP_GAM2_INTEGRATED 161

// See equations for free-free in gasphysics document, or in Rybicki and Lightman equations 5.14a
// and 5.18a
const double e6 = Constant::ESQUARE * Constant::ESQUARE * Constant::ESQUARE;
const double c3 = Constant::LIGHT * Constant::LIGHT * Constant::LIGHT;
const double sqrt_2piOver3m = sqrt(2. * Constant::PI / 3. / Constant::ELECTRONMASS);

// 32pi e^6 / 3mc^3 * sqrt(2pi / 3m)
const double gamma_nu_constantFactor =
                32 * Constant::PI * e6 / 3. / Constant::ELECTRONMASS / c3 * sqrt_2piOver3m;
// 4 e^6 / 3mhc * sqrt(2pi / 3m)
const double opCoef_nu_constantFactor = 4 * e6 / 3. / Constant::ELECTRONMASS / Constant::PLANCK /
                                        Constant::LIGHT * sqrt_2piOver3m;

// Constant factor from 1998-Sutherland equation 18
const double fk = 1.42554e-27;

FreeFree::FreeFree(const Array& frequencyv) : _frequencyv(frequencyv)
{
	readFullData(REPOROOT "/dat/gauntff_merged_Z01.dat");
	readIntegratedData(REPOROOT "/dat/gauntff_freqint_Z01.dat");
	DEBUG("Constructed FreeFree" << endl);
}

void FreeFree::readFullData(const string& file)
{
	// Translated to c++ from interpolate3.c that came with the 2014 van Hoof paper (MNRAS 444
	// 420)
	ifstream input(IOTools::ifstreamFile(file));

	// buffer
	string line;

	/* skip copyright statement */
	while (getline(input, line))
	{
		if (line.at(0) != '#')
			break;
	}

	/* read magic number */
	const long gaunt_magic = 20140210L;
	const long gaunt_magic2 = 20140510L;
	long magic;
	istringstream(line) >> magic;
	if (magic != gaunt_magic && magic != gaunt_magic2)
	{
		cerr << "read_table() found wrong magic number in file %s.\n" << endl;
		Error::runtime("read_table() found wrong magic number in file %s.\n");
	}

	/* read dimensions of the table */
	size_t np_gam2, np_u;
	getline(input, line);
	istringstream(line) >> np_gam2 >> np_u;
	assert(np_gam2 == NP_GAM2 && np_u == NP_U);

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
	if (magic == gaunt_magic2)
		getline(input, line);

	/* next lines are comments */
	while (getline(input, line))
	{
		if (line.at(0) != '#')
			break;
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

#ifdef DEBUG_CONTINUUM_DATA
	ofstream out;
	out.open("gauntff.dat");
	for (size_t ipu = 0; ipu < np_u; ipu++)
	{
		for (size_t ipg2 = 0; ipg2 < np_gam2; ++ipg2)
		{
			out << _fileGauntFactorvv(ipu, ipg2) << '\t';
		}
		out << endl;
	}
	out.close();
#endif
	DEBUG("Successfully read gauntff.dat" << endl);
}

void FreeFree::readIntegratedData(const string& file)
{
	// Translated to c++ from interpolate3.c that came with the 2014 van Hoof paper (MNRAS 444
	// 420)
	ifstream input(IOTools::ifstreamFile(file));

	// buffer
	string line;

	/* skip copyright statement */
	while (getline(input, line))
	{
		if (line.at(0) != '#')
			break;
	}

	/* read magic number */
	const long gaunt_magic = 20141008L;
	long magic;
	istringstream(line) >> magic;
	if (magic != gaunt_magic)
	{
		cerr << "read_table() found wrong magic number in file %s.\n" << endl;
		Error::runtime("read_table() found wrong magic number in file %s.\n");
	}

	/* read dimensions of the table */
	size_t np_gam2;
	getline(input, line);
	istringstream(line) >> np_gam2;
	assert(np_gam2 == NP_GAM2_INTEGRATED);

	/* read start value for log(gamma^2) */
	getline(input, line);
	istringstream(line) >> _loggamma2Min_integrated;

	/* read step size in dex */
	getline(input, line);
	istringstream(line) >> _logStep_integrated;
	_loggamma2Max_integrated = _loggamma2Min_integrated +
	                           static_cast<double>(np_gam2 - 1) * _logStep_integrated;

	/* read atomic number */
	getline(input, line);

	/* next lines are comments */
	while (getline(input, line))
	{
		if (line.at(0) != '#')
			break;
	}

	_fileGauntFactorv_integrated.resize(np_gam2);
	_loggamma2_integrated.resize(np_gam2);
	for (size_t ipg2 = 0; ipg2 < np_gam2; ++ipg2)
	{
		istringstream(line) >> _loggamma2_integrated[ipg2] >>
		                _fileGauntFactorv_integrated[ipg2];
		_fileGauntFactorv_integrated[ipg2] = log(_fileGauntFactorv_integrated[ipg2]);
		getline(input, line);
	}

#ifdef DEBUG_CONTINUUM_DATA
	ofstream out;
	out.open("integratedgauntff.dat");
	for (double logg2 = -5.9; logg2 < 9.9; logg2 += .1)
	{
		out << logg2 << '\t' << integratedGauntFactor(logg2) << '\n';
	}
	out.close();
#endif
}

double FreeFree::gauntFactor(double logu, double logg2) const
{
	// Throw an error if out of range for now. Maybe allow extrapolation later.
	Error::rangeCheck("log(gamma^2)", logg2, _loggamma2Min, _loggamma2Max);
	Error::rangeCheck("log(u)", logu, _loguMin, _loguMax);

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

	double gff = TemplatedUtils::interpolateRectangular<double>(
	                logg2, logu, xLeft, xRight, yLow, yUp, _fileGauntFactorvv(iLower, iLeft),
	                _fileGauntFactorvv(iLower, iRight), _fileGauntFactorvv(iUpper, iLeft),
	                _fileGauntFactorvv(iUpper, iRight));
	return exp(gff);
}

double FreeFree::integratedGauntFactor(double logg2) const
{
	// Throw an error if out of range for now. Maybe allow extrapolation later.
	Error::rangeCheck("log(gamma^2)", logg2, _loggamma2Min_integrated,
	                  _loggamma2Max_integrated);

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
	double loggff_integrated_ip =
	                TemplatedUtils::interpolateLinear(logg2, xLeft, xRight, fLeft, fRight);
	return exp(loggff_integrated_ip);
}

void FreeFree::addEmissionCoefficientv(double T, Array& gamma_nuv) const
{
	// gamma is fixed for a given temperature
	double kT = Constant::BOLTZMAN * T;
	double sqrtkT = sqrt(kT);
	double logg2 = log10(Constant::RYDBERG / kT);

#ifdef DEBUG_CONTINUUM_DATA
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/gammanuff.dat");
#endif
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double u = Constant::PLANCK * _frequencyv[iFreq] / kT;
		double logu = log10(u);
		double gammaNu = gamma_nu_constantFactor / sqrtkT * exp(-u) *
		                 gauntFactor(logu, logg2);
#ifdef DEBUG_CONTINUUM_DATA
		out << _frequencyv[iFreq] << "\t" << gammaNu / 1.e-40 << endl;
#endif
		gamma_nuv[iFreq] += gammaNu;
	}
#ifdef DEBUG_CONTINUUM_DATA
	out.close();
#endif
}

void FreeFree::addOpacityCoefficientv(double T, Array& opCoeffv) const
{
	double kT = Constant::BOLTZMAN * T;
	double sqrtkT = sqrt(kT);
	double loggamma2 = log10(Constant::RYDBERG / kT);
#ifdef DEBUG_CONTINUUM_DATA
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/ffopacitycoef.dat");
#endif
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		double nu = _frequencyv[iFreq];
		double u = Constant::PLANCK * nu / kT;
		double logu = log10(u);
		// C / nu^3 (1 - exp(-u)) gff(u, gamma^2)
		double opCoeffNu = opCoef_nu_constantFactor / sqrtkT / nu / nu / nu * -expm1(-u) *
		                   gauntFactor(logu, loggamma2);
#ifdef DEBUG_CONTINUUM_DATA
		out << _frequencyv[iFreq] << "\t" << opCoeffNu << endl;
#endif
		opCoeffv[iFreq] += opCoeffNu;
	}
#ifdef DEBUG_CONTINUUM_DATA
	out.close();
#endif
}

double FreeFree::heating(double np_ne, double T, const Array& specificIntensityv) const
{
	Array freefreeOpCoefv(_frequencyv.size());
	addOpacityCoefficientv(T, freefreeOpCoefv);
	Array intensityOpacityv(specificIntensityv * np_ne * freefreeOpCoefv);
	return Constant::FPI * TemplatedUtils::integrate<double>(_frequencyv, intensityOpacityv);
}

double FreeFree::cooling(double np_ne, double T) const
{
#ifdef DEBUG_CONTINUUM_DATA
	Array gamma_nuv(_frequencyv.size());
	addEmissionCoefficientv(T, gamma_nuv);

	// emissivity = ne np / 4pi * gamma
	// total emission = 4pi integral(emissivity) = ne np integral(gamma)
	double manually_integrated =
	                np_ne * TemplatedUtils::integrate<double>(_frequencyv, gamma_nuv);
#endif
	// Interpolate and use the frequency integrated gaunt factor
	double logg2 = log10(Constant::RYDBERG / Constant::BOLTZMAN / T);
	double from_data = np_ne * fk * sqrt(T) * integratedGauntFactor(logg2);
#ifdef DEBUG_CONTINUUM_DATA
	printf("Manually integrated: %e; from data: %e; ratio: %e\n", manually_integrated,
	       from_data, from_data / manually_integrated);
#endif
	return from_data;
}
