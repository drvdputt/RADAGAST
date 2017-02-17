#ifndef _TESTING_H_
#define _TESTING_H_
#include <vector>

namespace Testing
{
	std::vector<double> generateWavelengthGrid(size_t Nlambda, double lambdaMin, double lambdaMax);
	void refineWavelengthGrid(std::vector<double>& grid, size_t nPerLine, double spacingPower, std::vector<double> lineWavev, std::vector<double> lineWidthv);
	std::vector<double> generateISRF(const std::vector<double>& wavelength, double Tc, double G0);

	void testTwoLevel();
	void testGasSpecies();
}

#endif
