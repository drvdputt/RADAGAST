#include "Testing.h"
#include "Constants.h"
#include "IonizationBalance.h"
#include "NumUtils.h"
#include "PhotoelectricHeating.h"
#include "TemplatedUtils.h"
#include "TwoLevel.h"

#include "GasInterface.h"
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
                                  Array lineFreqv, Array freqWidthv)
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
		I_nu[frequencyv.size() - iWav - 1] = I_lambda[iWav] * wavelengthv[iWav] *
		                                     wavelengthv[iWav] / Constant::LIGHT;

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
	vector<double> wavelengthUV(wavelengthv.begin() + startLambdaUV,
	                            wavelengthv.begin() + endLambdaUV);
	vector<double> isrfUV(I_lambda.begin() + startLambdaUV, I_lambda.begin() + endLambdaUV);

	// Integrate over the UV only
	double UVdensity = Constant::FPI / Constant::LIGHT *
	                   NumUtils::integrate<double>(wavelengthUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING;

	// Rescale to _G0
	I_nu *= G0 / currentG0;

	vector<double> frequencyUV(end(frequencyv) - endLambdaUV - 1,
	                           end(frequencyv) - startLambdaUV - 1);
	vector<double> isrfUVbis(end(I_nu) - endLambdaUV - 1, end(I_nu) - startLambdaUV - 1);

	// Integrate over the UV only
	double UVdensitybis = Constant::FPI / Constant::LIGHT *
	                      NumUtils::integrate<double>(frequencyUV, isrfUVbis);
	cout << "Normalized spectrum uv = " << UVdensitybis << " ("
	     << UVdensitybis / Constant::HABING << " habing)" << endl;

	// Write out the ISRF
	std::ofstream out;
	out.open("isrf.txt");
	for (size_t b = 0; b < I_nu.size(); b++)
		out << frequencyv[b] << '\t' << I_nu[b] << '\n';
	out.close();
	out.open("isrfUV.txt");
	for (size_t b = 0; b < isrfUV.size(); b++)
		out << frequencyUV[b] << '\t' << isrfUVbis[b] << '\n';
	out.close();
	return I_nu;
}

Array Testing::freqToWavSpecificIntensity(const vector<double>& frequencyv,
                                          const Array& specificIntensity_nu)
{
	Array I_lambda(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
	{
		I_lambda[I_lambda.size() - iFreq - 1] = specificIntensity_nu[iFreq] *
		                                        frequencyv[iFreq] * frequencyv[iFreq] /
		                                        Constant::LIGHT;
	}
	return I_lambda;
}

void Testing::testIonizationStuff()
{
	ofstream out;

	double A0 = 6.30e-18;
	out.open("ionizationCrosSection.dat");
	for (double freq = 0; freq < 6.e16; freq += 1e14)
	{
		double sigma = Ionization::crossSection(freq);
		double eps = sqrt(freq / Ionization::THRESHOLD - 1);
		double sigmaTheoretical =
		                freq > Ionization::THRESHOLD
		                                ? A0 * pow(Ionization::THRESHOLD / freq, 4.) *
		                                                  exp(4 - 4 * atan(eps) / eps) /
		                                                  (1 - exp(-2 * Constant::PI / eps))
		                                : 0;
		out << freq << "\t" << sigma << "\t" << sigmaTheoretical << "\t"
		    << (sigma - sigmaTheoretical) / sigma << endl;
	}
	out.close();

	out.open("recombinationCooling.dat");
	out << "#kT / eV \t alpha" << endl;
	// by choosing these parameters, panel 2 of figure 9 from 2017-Mao should be reproduced
	double nH = 1;
	double f = 1;
	vector<double> kT_eVv = generateGeometricGridv(300, 1e-3, 1e3);
	for (double kT_eV : kT_eVv)
	{
		double T = kT_eV / Constant::BOLTZMAN / Constant::ERG_EV;
		double cool = Ionization::cooling(nH, f, T) / Constant::RYDBERG;
		out << kT_eV << "\t" << cool << endl;
	}
	out.close();

	// out.open("collisionalIonization")
}

void Testing::testGasInterfaceImpl()
{
	double Tc = 20000;
	double G0 = 1e0;
	double n = 1e1;
	double expectedTemperature = 1000;

	vector<double> tempFrequencyv =
	                generateGeometricGridv(10000, Constant::LIGHT / (1e10 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.00001 * Constant::UM_CM));

	//	const double lineWindowFactor = 5;
	//	vector<double> lineFreqv = {Constant::LIGHT * 82258.9191133};
	//	vector<double> decayRatev = {6.2649e+08};
	//	double thermalFactor =
	//	                sqrt(Constant::BOLTZMAN * 500000 / Constant::HMASS_CGS) /
	// Constant::LIGHT; 	vector<double> lineWidthv;
	//	lineWidthv.reserve(lineFreqv.size());
	//	for (size_t l = 0; l < lineFreqv.size(); l++)
	//	{
	//		double freq = lineFreqv[l];
	//		cout << "thermal width " << freq * thermalFactor << " natural width "
	//		     << decayRatev[l] << endl;
	//		lineWidthv.push_back(lineWindowFactor * (freq * thermalFactor +
	// decayRatev[l]));
	//	}
	//	refineFrequencyGrid(tempFrequencyv, 101, 3., lineFreqv, lineWidthv);

	Array frequencyv(tempFrequencyv.data(), tempFrequencyv.size());
	GasInterface gi(frequencyv, true);
	frequencyv = gi.frequencyv();

	Array specificIntensityv = generateSpecificIntensityv(
	                vector<double>(begin(frequencyv), end(frequencyv)), Tc, G0);

	GasState gs;
	gi.updateGasState(gs, n, expectedTemperature, specificIntensityv);

	const Array& emv = gs._emissivityv;
	const Array& opv = gs._opacityv;
	const Array& scav = gs._scatteringOpacityv * gs._previousISRFv;

	cout << "Integrated emissivity " << TemplatedUtils::integrate<double>(frequencyv, emv)
	     << endl;

	ofstream out, wavfile;
	char tab = '\t';
	out.open("opticalProperties.dat");
	vector<std::string> colnames = {
	                "frequency",
	                "wavelength",
	                "intensity j_nu (erg s-1 cm-3 Hz-1 sr-1)",
	                "opacity alpha_nu (cm-1)",
	                "scattered (erg s-1 cm-3 Hz-1 sr-1)",
	};
	out << "#";
	int i = 0;
	for (const auto& s : colnames)
	{
		out << i << ":" << s << tab;
		i++;
	}
	out << endl;
	wavfile.open("wavelengths.dat");
	wavfile << "#wav (micron)" << tab << "freq (Hz)" << endl;
	for (size_t iFreq = 0; iFreq < emv.size(); iFreq++)
	{
		double freq = frequencyv[iFreq];
		double wav = Constant::LIGHT / freq * Constant::CM_UM;
		out.precision(9);
		double effective = emv[iFreq] - scav[iFreq];
		effective = effective > 0 ? effective : 0;
		out << scientific << freq << tab << wav << tab << emv[iFreq] << tab << opv[iFreq]
		    << tab << scav[iFreq] << tab << emv[iFreq] - scav[iFreq] << tab << endl;
		wavfile.precision(9);
		wavfile << wav << tab << freq << endl;
	}
	out.close();
	wavfile.close();

	// Print some line intensities relative to Hbeta

	double fHalpha = Constant::LIGHT / 656.453e-7;
	double fHbeta = Constant::LIGHT / 486.264e-7;
	double fHgamma = Constant::LIGHT / 434.165e-7;

	double fLya = Constant::LIGHT / 121.567e-7;

	double fPalpha = Constant::LIGHT / 1875.61e-7;
	double fPbeta = Constant::LIGHT / 1282.16e-7;

	double fBralpha = Constant::LIGHT / 4052.27e-7;

	function<double(double frequency)> evaluateSpectrum = [&](double f) {
		return TemplatedUtils::evaluateLinInterpf(f, frequencyv, emv);
	};

	double Hbeta = evaluateSpectrum(fHbeta);
	cout << "Halpha / Hbeta " << evaluateSpectrum(fHalpha) / Hbeta << endl;
	cout << "Hgamma / Hbeta " << evaluateSpectrum(fHgamma) / Hbeta << endl;

	cout << "Lyalpha / Hbeta " << evaluateSpectrum(fLya) / Hbeta << endl;

	cout << "Palpha / Hbeta " << evaluateSpectrum(fPalpha) / Hbeta << endl;
	cout << "Pbeta / Hbeta " << evaluateSpectrum(fPbeta) / Hbeta << endl;

	cout << "Bralpha / HBeta " << evaluateSpectrum(fBralpha) / Hbeta << endl;

	cout << "TestHydrogenCalculator done" << endl;

	//	cout << "----------------------------------" << endl;
	//	cout << "plotting heating curve..." << endl;
	//	gi.testHeatingCurve(n, specificIntensityv);
}

void Testing::testPhotoelectricHeating()
{
	PhotoelectricHeatingRecipe phr;
	double T = 1000;
	vector<double> G0values;
	if (T == 1000)
	{
		G0values = {2.45e-2, 2.45e-1, 2.45e0, 2.45e1, 2.45e2};
	}
	if (T == 100)
	{
		G0values = {.75e-1, .75e0, .75e1, .75e2, .75e3};
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
