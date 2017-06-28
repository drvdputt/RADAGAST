#include "HydrogenLevels.h"
#include "Constants.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <iostream>
#include <vector>

HydrogenLevels::HydrogenLevels(bool hardcoded)
                : NLevel(chooseDataProvider(hardcoded))
{
}

LevelDataProvider* HydrogenLevels::chooseDataProvider(bool hardcoded) const
{
	if (hardcoded) return new HydrogenHardcoded();
	else return new HydrogenFromFiles(5);
}

Array HydrogenLevels::emissivityv(const Solution& s) const
{
	return lineEmissivityv(s) + boundBoundContinuum(s);
}

Array HydrogenLevels::boundBoundContinuum(const Solution& s) const
{
	// TODO: make sure that the right indices are obtained
	int index2s = 1;
	Array result(frequencyv().size());
	// 1984-Nussbaumer
	double constFactor = Constant::PLANCK / Constant::FPI * s.nv(index2s);
	double nu0 = (ev(index2s) - ev(0)) / Constant::PLANCK;
	double C = 202.0; // s-1
	double alpha = .88;
	double beta = 1.53;
	double gam = .8;
	for (size_t iFreq = 0; frequencyv()[iFreq] < nu0; iFreq++)
	{
		double y = frequencyv()[iFreq] / nu0;
		double y1miny = y * (1 - y);
		double pow4y1miny_gam = pow(4 * y1miny, gam);
		double Py = C * (y1miny * (1 - pow4y1miny_gam) +
		                 alpha * pow(y1miny, beta) * pow4y1miny_gam);
		result[iFreq] = constFactor * y * Py;
	}
	return result;
}
