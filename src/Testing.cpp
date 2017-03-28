#include "Constants.h"
#include "HydrogenCalculator.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "PhotoelectricHeating.h"
#include "TemplatedUtils.h"
#include "Testing.h"
#include "TwoLevel.h"

#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>

using namespace std;

vector<double> Testing::generateGeometricGridv(size_t nPoints, double min, double max)
{
	vector<double> frequencyv(nPoints);
	double freqStepFactor = std::pow(max / min, 1. / (nPoints - 1));
	double freq = min;
	for (size_t n = 0; n < nPoints; n++)
	{
		frequencyv[n] = freq;
		freq *= freqStepFactor;
	}
	return frequencyv;
}

vector<double> Testing::freqToWavGrid(const vector<double>& frequencyv)
{
	vector<double> wavelengthv;
	wavelengthv.reserve(frequencyv.size());
	for (auto rit = frequencyv.rbegin(); rit != frequencyv.rend(); rit++)
		wavelengthv.push_back(Constant::LIGHT / *rit);
	return wavelengthv;
}

void Testing::refineFrequencyGrid(vector<double>& grid, size_t nPerLine, double spacingPower,
		vector<double> lineFreqv, vector<double> freqWidthv)
{
	// We want an odd nPerLine, so we can put 1 point in the center
	if (!(nPerLine % 2))
		nPerLine += 1;

	grid.reserve(grid.size() + nPerLine * lineFreqv.size());

	for (size_t i = 0; i < lineFreqv.size(); i++)
	{
		// Add a point at the center of the line, while keeping the vector sorted
		TemplatedUtils::sortedInsert<double, vector<double>>(lineFreqv[i], grid);

		// Add the rest of the points in a power law spaced way
		if (nPerLine > 1)
		{
			size_t nOneSide = (nPerLine - 1) / 2;
			double a = freqWidthv[i] / pow(nOneSide, spacingPower);
			for (size_t sidePoint = 1; sidePoint <= nOneSide; sidePoint++)
			{
				double distance = a * pow(sidePoint, spacingPower);

				// Left of center
				double freq = lineFreqv[i] - distance;
				TemplatedUtils::sortedInsert<double>(freq, grid);

				// Right of center
				freq = lineFreqv[i] + distance;
				TemplatedUtils::sortedInsert<double>(freq, grid);
			}
		}
	}
}

Array Testing::generateSpecificIntensityv(const vector<double>& frequencyv, double Tc, double G0)
{
	// A blackbody (in specific intensity per wavelength units)
	vector<double> wavelengthv = freqToWavGrid(frequencyv);
	vector<double> I_lambda = NumUtils::bbodyCGS<double>(wavelengthv, Tc);

	// Convert to per frequency units using I_nu = I_lambda * lambda * lambda / c
	Array I_nu(frequencyv.size());
	for (size_t iWav = 0; iWav < wavelengthv.size(); iWav++)
		I_nu[frequencyv.size() - iWav - 1] = I_lambda[iWav] * wavelengthv[iWav] * wavelengthv[iWav]
				/ Constant::LIGHT;

	// Cut out the UV part
	size_t i = 0;
	size_t startLambdaUV, endLambdaUV;
	while (wavelengthv[i] < 912 * Constant::ANG_CM && i < I_lambda.size())
		i++;
	startLambdaUV = i > 0 ? i - 1 : 0;
	while (wavelengthv[i] < 2400 * Constant::ANG_CM && i < I_lambda.size())
		i++;
	endLambdaUV = i + 1;
	cout << "UV goes from " << startLambdaUV << " to " << endLambdaUV << endl;
	vector<double> wavelengthUV(wavelengthv.begin() + startLambdaUV, wavelengthv.begin() + endLambdaUV);
	vector<double> isrfUV(I_lambda.begin() + startLambdaUV, I_lambda.begin() + endLambdaUV);

	// Integrate over the UV only
	double UVdensity = Constant::FPI / Constant::LIGHT
			* NumUtils::integrate<double>(wavelengthUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING;

	// Rescale to _G0
	I_nu *= G0 / currentG0;

	vector<double> frequencyUV(end(frequencyv) - endLambdaUV - 1, end(frequencyv) - startLambdaUV - 1);
	vector<double> isrfUVbis(end(I_nu) - endLambdaUV - 1, end(I_nu) - startLambdaUV - 1);

	// Integrate over the UV only
	double UVdensitybis = -Constant::FPI / Constant::LIGHT
			* NumUtils::integrate<double>(frequencyUV, isrfUVbis);
	cout << "Normalized spectrum uv = " << UVdensitybis << " (" << UVdensitybis / Constant::HABING
			<< " habing)" << endl;

	// Write out the ISRF
	std::ofstream out;
	out.open("/Users/drvdputt/GasModule/run/isrf.txt");
	for (size_t b = 0; b < I_nu.size(); b++)
		out << frequencyv[b] << '\t' << I_nu[b] << '\n';
	out.close();
	out.open("/Users/drvdputt/GasModule/run/isrfUV.txt");
	for (size_t b = 0; b < isrfUV.size(); b++)
		out << frequencyUV[b] << '\t' << isrfUVbis[b] << '\n';

	return I_nu;
}

Array Testing::freqToWavSpecificIntensity(const vector<double>& frequencyv, const Array& specificIntensity_nu)
{
	Array I_lambda(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
	{
		I_lambda[I_lambda.size() - iFreq - 1] = specificIntensity_nu[iFreq] * frequencyv[iFreq]
				* frequencyv[iFreq] / Constant::LIGHT;
	}
	return I_lambda;
}

void Testing::testIonizationCrossSection()
{
	ofstream out;
	out.open("/Users/drvdputt/GasModule/run/ionizationCrosSection.dat");
	for (double freq = 0; freq < 6.e16; freq += 1e14)
	{
		double sigma = Ionization::crossSection(freq);
		out << freq << "\t" << sigma << endl;
	}
	out.close();
}

void Testing::testHydrogenCalculator()
{
	double Tc = 30000;
	double G0 = 1e0;
	double n = 1e0;
	double expectedTemperature = 6000;

	vector<double> tempFrequencyv = generateGeometricGridv(200,
			Constant::LIGHT / (1000 * Constant::UM_CM),
			Constant::LIGHT / (0.01 * Constant::UM_CM));

	const double lineWindowFactor = 5;
	vector<double> lineFreqv =
	{ Constant::LIGHT * 82258.9191133 };
	vector<double> decayRatev =
	{ 6.2649e+08 };
	double thermalFactor = sqrt(Constant::BOLTZMAN * 500000 / Constant::HMASS_CGS) / Constant::LIGHT;
	vector<double> lineWidthv;
	lineWidthv.reserve(lineFreqv.size());
	for (size_t l = 0; l < lineFreqv.size(); l++)
	{
		double freq = lineFreqv[l];
		cout << "thermal width " << freq * thermalFactor << " natural width " << decayRatev[l]
				<< endl;
		lineWidthv.push_back(lineWindowFactor * (freq * thermalFactor + decayRatev[l]));
	}
	refineFrequencyGrid(tempFrequencyv, 101, 3., lineFreqv, lineWidthv);

	Array frequencyv(tempFrequencyv.data(), tempFrequencyv.size());
	Array specificIntensityv = generateSpecificIntensityv(tempFrequencyv, Tc, G0);

	HydrogenCalculator hc(frequencyv);
	hc.solveBalance(n, expectedTemperature, specificIntensityv);

	const Array& lumv = hc.emissivityv();
	const Array& opv = hc.opacityv();
	const Array& scav = hc.scatteredv();

	cout << "Integrated emissivity " << TemplatedUtils::integrate<double>(frequencyv, lumv) << endl;

	ofstream out, wavfile;
	char tab = '\t';
	out.open("/Users/drvdputt/GasModule/run/opticalProperties.dat");
	out << "# 0:frequency" << tab << "1:intensity j_nu (erg s-1 cm-3 Hz-1 sr-1)" << tab
			<< "2:opacity alpha_nu (cm-1)" << tab << "3:scattered (erg s-1 cm-3 Hz-1 sr-1)" << tab
			<< "4: int - sca" << endl;
	wavfile.open("/Users/drvdputt/GasModule/run/wavelengths.dat");
	for (size_t iFreq = 0; iFreq < lumv.size(); iFreq++)
	{
		double freq = frequencyv[iFreq];
		double wav = Constant::LIGHT / freq * Constant::CM_UM;
		out.precision(9);
		double effective = lumv[iFreq] - scav[iFreq];
		effective = effective > 0 ? effective : 0;
		out << scientific << freq << tab << lumv[iFreq] << tab << opv[iFreq] << tab << scav[iFreq]
				<< tab << lumv[iFreq] - scav[iFreq] << tab << effective << endl;
		wavfile.precision(9);
		wavfile << wav << endl;
	}
	out.close();
	wavfile.close();

	cout << "----------------------------------" << endl;
	cout << "plotting heating curve..." << endl;
	hc.testHeatingCurve();
}

void Testing::testPhotoelectricHeating()
{
	PhotoelectricHeatingRecipe phr;
	double T = 1000;
	vector<double> G0values;
	if (T == 1000)
	{
		G0values =
		{	2.45e-2, 2.45e-1, 2.45e0, 2.45e1, 2.45e2};
	}
	if (T == 100)
	{
		G0values =
		{	.75e-1, .75e0, .75e1, .75e2, .75e3};
	}

	phr.setGasTemperature(T);
	for (double G0 : G0values)
	{
		phr.setG0(G0);
		stringstream filename;
		filename << "/Users/drvdputt/GasModule/run/photoelectricHeatingG" << setprecision(4)
				<< scientific << G0 << ".dat";
		phr.heatingRateTest(filename.str());
	}
}
