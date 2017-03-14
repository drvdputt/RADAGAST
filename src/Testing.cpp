#include "Constants.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "Testing.h"
#include "TwoLevel.h"

#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "HydrogenCalculator.h"
#include "TemplatedUtils.h"
#include "PhotoelectricHeating.h"

using namespace std;

vector<double> Testing::generateFrequencyGrid(size_t nFreq, double minFreq, double maxFreq)
{
	vector<double> frequencyv(nFreq);
	double freqStepFactor = std::pow(maxFreq / minFreq, 1. / (nFreq - 1));
	double freq = minFreq;
	for (size_t n = 0; n < nFreq; n++)
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
		TemplatedUtils::sorted_insert<double>(grid, lineFreqv[i]);

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
				TemplatedUtils::sorted_insert<double>(grid, freq);

				// Right of center
				freq = lineFreqv[i] + distance;
				TemplatedUtils::sorted_insert<double>(grid, freq);
			}
		}
	}
}

std::vector<double> Testing::generateSpecificIntensity(const std::vector<double>& frequencyv, double Tc,
		double G0)
{
	// A blackbody (in specific intensity per wavelength units)
	vector<double> wavelengthv;
	wavelengthv.reserve(frequencyv.size());
	for (auto rit = frequencyv.rbegin(); rit != frequencyv.rend(); rit++)
		wavelengthv.push_back(Constant::LIGHT / *rit);

	vector<double> I_lambda = NumUtils::bbodyCGS<double>(wavelengthv, Tc);

	// Convert to per frequency units using I_nu = I_lambda * lambda * lambda / c
	vector<double> I_nu;
	I_nu.reserve(frequencyv.size());
	auto wavRit = wavelengthv.rbegin();
	for (auto IlambdaRit = I_lambda.rbegin(); IlambdaRit < I_lambda.rend(); IlambdaRit++, wavRit++)
		I_nu.push_back((*IlambdaRit) * (*wavRit) * (*wavRit) / Constant::LIGHT);

	// Cut out the UV part
	size_t i = 0;
	size_t startUV, endUV;
	while (wavelengthv[i] < 912 * Constant::ANG_CM && i < I_lambda.size())
		i++;
	startUV = i > 0 ? i - 1 : 0;
	while (wavelengthv[i] < 2400 * Constant::ANG_CM && i < I_lambda.size())
		i++;
	endUV = i + 1;
	cout << "UV goes from " << startUV << " to " << endUV << endl;
	vector<double> wavelengthUV(wavelengthv.begin() + startUV, wavelengthv.begin() + endUV);
	vector<double> isrfUV(I_lambda.begin() + startUV, I_lambda.begin() + endUV);

	// Integrate over the UV only
	double UVdensity = Constant::FPI / Constant::LIGHT
			* NumUtils::integrate<double>(wavelengthUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING;

	// Rescale to _G0
	for (double& d : I_nu)
		d *= G0 / currentG0;

	vector<double> frequencyUV(frequencyv.rbegin() + startUV, frequencyv.rbegin() + endUV);
	vector<double> isrfUVbis(I_nu.rbegin() + startUV, I_nu.rbegin() + endUV);

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

vector<double> Testing::freqToWavSpecificIntensity(const vector<double>& frequencyv,
		const vector<double>& specificIntensity_nu)
{
	vector<double> I_lambda;
	I_lambda.reserve(frequencyv.size());
	auto fRit = frequencyv.rbegin();
	auto InuRit = specificIntensity_nu.rbegin();
	while (InuRit != specificIntensity_nu.rend())
	{
		double freq = *fRit;
		I_lambda.push_back(*InuRit * freq * freq / Constant::LIGHT);
		InuRit++;
		fRit++;
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
	double expectedTemperature = 10000;

	vector<double> frequencyv = generateFrequencyGrid(2000, Constant::LIGHT / (200 * Constant::UM_CM),
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
	refineFrequencyGrid(frequencyv, 1001, 3., lineFreqv, lineWidthv);

	vector<double> specificIntensity = generateSpecificIntensity(frequencyv, Tc, G0);

	HydrogenCalculator hc(frequencyv);
	hc.solveBalance(n, expectedTemperature, specificIntensity);

	const vector<double>& lumv = hc.emissivityv();
	const vector<double>& opv = hc.opacityv();
	const vector<double>& scav = hc.scatteredv();

	cout << "Integrated emissivity " << NumUtils::integrate<double>(frequencyv, lumv) << endl;

	ofstream out;
	char tab = '\t';
	out.open("/Users/drvdputt/GasModule/run/opticalProperties.dat");
	for (size_t iFreq = 0; iFreq < lumv.size(); iFreq++)
	{
		double freq = frequencyv[iFreq];
		//double wav = Constant::LIGHT / freq;
		out.precision(9);
		out << scientific << freq << tab << lumv[iFreq] << tab << opv[iFreq] << tab << scav[iFreq]
				<< tab << lumv[iFreq] - scav[iFreq] << endl;
	}
	out.close();

//	cout << "----------------------------------" << endl;
//	cout << "plotting heating curve..." << endl;
//	hc.testHeatingCurve();
}

void Testing::testPhotoelectricHeating()
{
	PhotoelectricHeatingRecipe phr;
	double T = 1000;
	vector<double> G0values;
	if (T == 1000)
	{
		G0values =
		{ 2.45e-2, 2.45e-1, 2.45e0, 2.45e1, 2.45e2 };
	}
	if (T == 100)
	{
		G0values =
		{ .75e-1, .75e0, .75e1, .75e2, .75e3};
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
