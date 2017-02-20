#include "Testing.h"
#include "TwoLevel.h"
#include "GasSpecies.h"

#include <Constants.h>
#include <NumUtils.h>
#include <iostream>
#include <fstream>
#include <ios>

using namespace std;

std::vector<double> Testing::generateWavelengthGrid(size_t Nlambda, double lambdaMin, double lambdaMax)
{
	vector<double> wavelength(Nlambda);
	double lambdaStepFactor = std::pow(lambdaMax / lambdaMin, 1. / Nlambda);
	float lambda = lambdaMin;
	for (size_t n = 0; n < Nlambda; n++)
	{
		wavelength[n] = lambda;
		lambda *= lambdaStepFactor;
	}
	return wavelength;
}

// Help function to refine the wavelength grid
namespace {
	template <typename T>
	void sorted_insert(vector<T>& vec, T elem)
	{
		vec.insert(upper_bound(vec.begin(), vec.end(), elem), elem);
	}
}

void Testing::refineWavelengthGrid(vector<double>& grid, size_t nPerLine, double spacingPower, vector<double> lineWavev, vector<double> lineWidthv)
{
	// We want an odd nPerLine, so we can put 1 point in the center
	if (!(nPerLine % 2)) nPerLine += 1;

	grid.reserve(grid.size() + nPerLine * lineWavev.size());

	for (size_t i = 0; i < lineWavev.size(); i++)
	{
		// Add a point at the center of the line, while keeping the vector sorted
		sorted_insert<double>(grid, lineWavev[i]);

		// Add the rest of the points in a power law spaced way
		if (nPerLine > 1)
		{
			cout << "Putting extra grid points at wavelengths ";
			size_t nOneSide = (nPerLine - 1) / 2;
			double a = lineWidthv[i] / pow(nOneSide, spacingPower);
			for (size_t sidePoint = 1; sidePoint <= nOneSide; sidePoint++)
			{
				double distance = a * pow(sidePoint, spacingPower);

				// Left of center
				double wav = lineWavev[i] - distance;
				cout << wav * Constant::CM_UM << " ";
				sorted_insert<double>(grid, wav);

				// Right of center
				wav = lineWavev[i] + distance;
				cout << wav * Constant::CM_UM << " ";
				sorted_insert<double>(grid, wav);
			}
			cout << endl;
		}
	}
}

std::vector<double> Testing::generateISRF(const std::vector<double>& wavelength, double Tc, double G0)
{
	// A blackbody
	vector<double> isrf = NumUtils::bbodyCGS<double>(wavelength, Tc);

	// Convert to energy density
	for (double& d : isrf) d *= Constant::FPI / Constant::LIGHT;

	// Cut out the UV part
	size_t i = 0;
	size_t startUV, endUV;
	while (wavelength[i] < 912 * Constant::ANG_CM && i < isrf.size()) i++;
	startUV = i > 0 ? i - 1 : 0;
	while (wavelength[i] < 2400 * Constant::ANG_CM && i < isrf.size()) i++;
	endUV = i + 1;

	vector<double> wavelengthUV(wavelength.begin() + startUV, wavelength.begin() + endUV);
	vector<double> isrfUV(isrf.begin() + startUV, isrf.begin() + endUV);
	cout << "UV goes from " << startUV << " (" << wavelengthUV[0] * Constant::CM_UM << ")"
		 << " to " << endUV << "("<< wavelengthUV[wavelengthUV.size() - 1] * Constant::CM_UM << ")" << endl;

	// Integrate over the UV only
	double UVdensity = NumUtils::integrate<double>(wavelengthUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING;

	// Rescale to _G0
	for (double& d : isrf) d *= G0 / currentG0;

	vector<double> wavelengthUVbis(wavelength.begin() + startUV, wavelength.begin() + endUV);
	vector<double> isrfUVbis(isrf.begin() + startUV, isrf.begin() + endUV);

	// Integrate over the UV only
	double UVIntensitybis = NumUtils::integrate<double>(wavelengthUVbis, isrfUVbis);
	cout << "Normalized spectrum uv = "<< UVIntensitybis << " (" <<  UVIntensitybis / Constant::HABING << " habing)"<< endl;

	// Write out the ISRF
	std::ofstream out;
	out.open("/Users/drvdputt/Testing/isrf.txt");
	for (size_t b = 0; b < isrf.size(); b++)
		out << wavelength[b] << '\t' << isrf[b] << '\n';
	out.close();
	out.open("/Users/drvdputt/Testing/isrfUV.txt");
	for (size_t b = 0; b < isrfUV.size(); b++)
		out << wavelengthUV[b] << '\t'<< isrfUV[b] << '\n';

	return isrf;
}

void Testing::testTwoLevel()
{
	double Tc = 6000;
	double G0 = 2000;

	const std::vector<double>& wavelength = generateWavelengthGrid(200, 157.6 * Constant::UM_CM, 157.9 * Constant::UM_CM);
	const std::vector<double>& isrf = generateISRF(wavelength, Tc, G0);

	double n = 25;
	TwoLevel tl(wavelength);

	tl.doLevels(n, n, 50000, isrf, 0);

	double lum = tl.bolometricEmission(1, 0);
	const std::vector<double>& lumv = tl.calculateEmission();
	const std::vector<double>& opv = tl.calculateOpacity();

	cout << "Total line emissivity " << lum << endl;
	cout << "Integrated emissivity " << NumUtils::integrate<double>(wavelength, lumv) << endl;

	ofstream em_out, op_out;
	em_out.open("/Users/drvdputt/GasModule/bin/emission.dat");
	op_out.open("/Users/drvdputt/GasModule/bin/opacity.dat");
	for(size_t w = 0; w < lumv.size(); w++)
	{
		em_out << wavelength[w] << '\t'<< lumv[w] << endl;
		op_out << wavelength[w] << '\t'<< opv[w] << endl;
	}
	em_out.close();
	op_out.close();
}

void Testing::testGasSpecies()
{
	double Tc = 10000;
	double G0 = 0.1;
	double n = 10.;
	double expectedTemperature = 1000;

	vector<double> wavelength = generateWavelengthGrid(200, 0.01 * Constant::UM_CM, 200 * Constant::UM_CM);
	vector<double> lineWavev = {157.740709 * Constant::UM_CM};
	double thermalFactor = sqrt(Constant::BOLTZMAN * expectedTemperature / Constant::HMASS_CGS / Constant::LIGHT/Constant::LIGHT);
	vector<double> lineWidthv;
	lineWidthv.reserve(lineWavev.size());
	for (double wav : lineWavev)
		lineWidthv.push_back(4.5 * wav * thermalFactor);
	refineWavelengthGrid(wavelength, 12, 2, lineWavev, lineWidthv);

	GasSpecies gs(wavelength);

	vector<double> isrf = generateISRF(wavelength, Tc, G0);

	gs.solveBalance(n, isrf);

	const vector<double>& lumv = gs.emissivity();
	const vector<double>& opv = gs.opacity();

	cout << "Integrated emissivity " << NumUtils::integrate<double>(wavelength, lumv) << endl;

	ofstream em_out, op_out;
	em_out.open("/Users/drvdputt/GasModule/bin/emission.dat");
	op_out.open("/Users/drvdputt/GasModule/bin/opacity.dat");
	for(size_t w = 0; w < lumv.size(); w++)
	{
		em_out.precision(9);
		em_out << scientific << wavelength[w] << '\t'<< lumv[w] << endl;
		op_out.precision(9);
		op_out << scientific << wavelength[w] << '\t'<< opv[w] << endl;
	}
	em_out.close();
	op_out.close();
}
