#include "FreeFree.h"
#include "Constants.h"
#include "flags.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

#define NP_GAM2 81
#define NP_U 146

FreeFree::FreeFree()
{
	readData("/Users/drvdputt/GasModule/git/dat/gauntff_merged_Z01.dat");
}

void FreeFree::readData(string file)
{
	// Translated to c++ from interpolate3.c that came with the 2014 van Hoof paper (MNRAS 444 420)
	ifstream input;
	input.open(file);
	if (!input)
		throw std::runtime_error("File " + file + "not found.");

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
		throw std::runtime_error("read_table() found wrong magic number in file %s.\n");
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

	_loggamma2Max = _loggamma2Min + (double) (np_gam2 - 1) * _logStep;
	_loguMax = _loguMin + (double) (np_u - 1) * _logStep;

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

#ifdef PRINT_CONTINUUM_DATA
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/gauntff.dat");
	for (size_t ipu=0; ipu < np_u; ipu++)
	{
		for (size_t ipg2=0; ipg2 < np_gam2; ++ipg2)
		{
			out << _fileGauntFactorvv(ipu, ipg2) << '\t';
		}
		out << endl;
	}
	out.close();
#endif
}

double FreeFree::gauntFactor(double temperature, double frequency) const
{
	// Since this data is not given for frequencies, but for the values u and gamma^2 instead,
	// we cannot simply precalculate one of the interpolations like we did for the free-bound.
	// (both parameters are temperature dependent)
	double kT = Constant::BOLTZMAN * temperature;
	double logg2 = log(Constant::RYDBERG / kT);
	double u = log(Constant::PLANCK * frequency / kT);

	// Interpolate bi-linearly between the data points
	return 0;
}
