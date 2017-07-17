#include "Testing.h"
#include "Constants.h"
#include "Error.h"
#include "FreeBound.h"
#include "GasInterface.h"
#include "GasInterfaceImpl.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "NLevel.h"
#include "NumUtils.h"
#include "PhotoelectricHeating.h"
#include "TemplatedUtils.h"

#include <Eigen/Dense>

using namespace std;

Array Testing::improveFrequencyGrid(const NLevel& boundBound, const FreeBound& freeBound,
                                    const Array& oldPoints)
{
	// Add extra points for the lines
	int numLines;
	Array lineFreqv, lineWidthv;
	boundBound.lineInfo(numLines, lineFreqv, lineWidthv);

	double lineWindowFactor = 1.;
	double thermalFactor =
	                sqrt(Constant::BOLTZMAN * 500000 / Constant::HMASS_CGS) / Constant::LIGHT;
	lineWidthv = lineWindowFactor * (lineWidthv + lineFreqv * thermalFactor);

	vector<double> gridVector(begin(oldPoints), end(oldPoints));
	Testing::refineFrequencyGrid(gridVector, 13, 2.5, lineFreqv, lineWidthv);

	// And for the jumps in the bound-bound spectrum
	const Array& thresholdv = freeBound.thresholdv();
	// Don't bother with the last jump, because that's the end of the data
	Array jumpFreqv(&thresholdv[0], thresholdv.size() - 2);
	Testing::refineFrequencyGrid(gridVector, 3, 1., jumpFreqv, 1e-6 * jumpFreqv);

	// And for the ionization threshold
	Array ionThr({Ionization::THRESHOLD});
	Testing::refineFrequencyGrid(gridVector, 3, 1., ionThr, 1e-6 * ionThr);

	return Array(gridVector.data(), gridVector.size());
}

vector<double> Testing::generateGeometricGridv(size_t nPoints, double min, double max)
{
	vector<double> frequencyv(nPoints);
	double freqStepFactor = pow(max / min, 1. / (nPoints - 1));
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

		// Skip the wing points of line if it lies near the core of the previous one, to prevent too much
		// points from bunching up.
		if (i > 0 && (lineFreqv[i] - lineFreqv[i - 1]) / freqWidthv[i - 1] < 0.01)
			continue;


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
	ofstream out = IOTools::ofstreamFile("isrf.txt");
	for (size_t b = 0; b < I_nu.size(); b++)
		out << frequencyv[b] << '\t' << I_nu[b] << '\n';
	out.close();
	out = IOTools::ofstreamFile("testing/isrfUV.txt");
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
	out = IOTools::ofstreamFile("ionization/crossSection.dat");
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

	out = IOTools::ofstreamFile("ionization/recombinationCooling.dat");
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
}

void Testing::runGasInterfaceImpl(const GasInterface& gi, const std::string& outputPath, double Tc,
                                  double G0, double n, double expectedTemperature)
{
	const Array& frequencyv = gi.frequencyv();

	Array specificIntensityv = generateSpecificIntensityv(
	                vector<double>(begin(frequencyv), end(frequencyv)), Tc, G0);

	GasState gs;
	gi.updateGasState(gs, n, expectedTemperature, specificIntensityv);

	const Array& emv = gs._emissivityv;
	const Array& opv = gs._opacityv;
	const Array& scav = gs._scatteringOpacityv * gs._previousISRFv;

	cout << "Integrated emissivity " << TemplatedUtils::integrate<double>(frequencyv, emv)
	     << endl;

	char tab = '\t';
	ofstream out = IOTools::ofstreamFile(outputPath + "opticalProperties.dat");
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
	ofstream wavfile = IOTools::ofstreamFile(outputPath + "wavelengths.dat");
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
		return TemplatedUtils::evaluateLinInterpf<double>(f, frequencyv, emv);
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
	//	plotHeatingCurve(*gi.pimpl(), outputPath, specificIntensityv, n);
}

void Testing::testPhotoelectricHeating()
{
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

	PhotoelectricHeatingRecipe phr;
	phr.setGasTemperature(T);

	phr.yieldFunctionTest();
	for (double G0 : G0values)
	{
		phr.chargeBalanceTest(G0);
		phr.heatingRateTest(G0);
	}
}

void Testing::testACollapse()
{
	HydrogenFromFiles hff(0);
	Eigen::MatrixXd avv = hff.avv().array();
	cout << hff.avv() << endl;
	cout << "Compare this with values directly from NIST below:" << endl;
	Eigen::MatrixXd nistA(5, 5);
	// clang-format off
	nistA << 0, 0, 0, 0, 0,
			4.6986e+08, 0, 0, 0, 0,
			5.5751e+07, 4.4101e+07, 0, 0, 0,
			1.2785e+07, 8.4193e+06, 8.9860e+06, 0, 0,
			4.1250e+06, 2.5304e+06, 2.2008e+06, 2.6993e+06, 0;
	// clang-format on
	cout << nistA << endl;
	cout << "The element-wise relative difference is " << endl;
	Eigen::MatrixXd relDiff = (avv - nistA).array() / nistA.array();
	// Take out all the nan's
	for (int i = 0; i < relDiff.size(); i++)
	{
		double* pointer = relDiff.data() + i;
		double value = *pointer;
		*pointer = isnan(value) ? 0 : value;
	}
	cout << relDiff << endl;
	// For automated testing, assert the following (maximum allowed error is just slightly above
	// the value outputted above)
	assert((relDiff.cwiseAbs().array() < 0.002).all());
}

void Testing::testPS64Collisions()
{
	const double T = 10000;
	const double ne = 1e4;
	const double np = 1e4;

	HydrogenFromFiles hff(5);
	Eigen::MatrixXd avv = hff.avv();
	Eigen::MatrixXd cvv = hff.cvv(T, ne, np);

	/* Calculate and write out (q_n(l-1) + q_n(l+1)) / A_nl, where A_nl is the total downwards
	   rate from level nl. */
	Eigen::VectorXd anlv = avv.rowwise().sum();
	for (int n : std::array<int, 2>{4, 5})
	{
		ofstream out = IOTools::ofstreamFile("ps64/t" + to_string(n) + "l_q" +
		                                     to_string(n) + "l.dat");
		for (int li = 0; li < n; li++)
		{
			int nliIndex = hff.indexOutput(n, li);
			// Decay rate to all other levels
			double anl = anlv(nliIndex);
			double qnl = 0;

			// Get the total decay rate due to l-changing collisions
			for (int lf = 0; lf < n; lf++)
			{
				int nlfIndex = hff.indexOutput(n, lf);
				qnl += cvv(nliIndex, nlfIndex);

				// l-changing collision rates should be zero except for changes by 1
				int deltal = li - lf;
				if (cvv(nliIndex, nlfIndex) > 0 && abs(deltal) != 1)
				{
					stringstream ss;
					ss << "The l-changing coefficient from " << n << "," << li
					   << " to " << n << "," << lf << " should be zero.";
					Error::runtime(ss.str());
				}
			}
			out << li << "\t" << qnl / anl << "\t" << qnl << "\t" << anl << endl;
		}
		out.close();
	}
}

void Testing::compareFromFilesvsHardCoded()
{
	HydrogenHardcoded hhc;
	HydrogenFromFiles hff(2);

	auto hc_vs_ff = [&](auto hc_thing, auto ff_thing) {
		cout << "Hardcoded" << endl;
		cout << hc_thing << endl;
		cout << "From files" << endl;
		cout << ff_thing << endl;
	};

	assert(hhc.numLv() == hff.numLv());

	cout << "Energy levels:" << endl;
	Eigen::VectorXd evhc = hhc.ev();
	Eigen::VectorXd evff = hff.ev();
	hc_vs_ff(evhc, evff);

	assert(hhc.gv() == hff.gv());

	cout << "A coefficients:" << endl;
	Eigen::MatrixXd avvhc = hhc.avv();
	Eigen::MatrixXd avvff = hff.avv();
	hc_vs_ff(avvhc, avvff);

	cout << "Extra A:" << endl;
	Eigen::MatrixXd eavvhc = hhc.extraAvv();
	Eigen::MatrixXd eavvff = hff.extraAvv();
	hc_vs_ff(eavvhc, eavvff);

	double T = 1e4;
	double ne = 1e4;
	double np = 1e4;

	cout << "Collisions:" << endl;
	Eigen::MatrixXd cvvhc = hhc.cvv(T, ne, np);
	Eigen::MatrixXd cvvff = hff.cvv(T, ne, np);
	hc_vs_ff(cvvhc, cvvff);

	cout << "Recombinations:" << endl;
	Eigen::VectorXd alphavhc = hhc.sourcev(T, ne, np) / ne / np;
	Eigen::VectorXd alphavff = hff.sourcev(T, ne, np) / ne / np;
	hc_vs_ff(alphavhc, alphavff);
}

void Testing::runFromFilesvsHardCoded()
{
	vector<double> tempFrequencyv =
	                generateGeometricGridv(1000, Constant::LIGHT / (1e10 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.00001 * Constant::UM_CM));
	Array unrefined(tempFrequencyv.data(), tempFrequencyv.size());

	// Hey, at least we'll get a decent frequency grid out of this hack
	HydrogenLevels hl(make_shared<HydrogenFromFiles>(5), unrefined);
	FreeBound fb(unrefined);
	Array frequencyv = improveFrequencyGrid(hl, fb, unrefined);

	GasInterface gihhc(frequencyv, "hhc");
	runGasInterfaceImpl(gihhc, "hardcoded/");

	GasInterface gihff(frequencyv, "hff2");
	runGasInterfaceImpl(gihff, "fromfiles/");
}

void Testing::runFullModel()
{
	vector<double> tempFrequencyv =
	                generateGeometricGridv(200, Constant::LIGHT / (1e4 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.005 * Constant::UM_CM));
	Array unrefined(tempFrequencyv.data(), tempFrequencyv.size());

	HydrogenLevels hl(make_shared<HydrogenFromFiles>(), unrefined);
	FreeBound fb(unrefined);
	Array frequencyv = improveFrequencyGrid(hl, fb, unrefined);

	GasInterface gihffFull(frequencyv, "hff");
	runGasInterfaceImpl(gihffFull, "");
}

void Testing::plotHeatingCurve(const GasInterfaceImpl& gi, const std::string& outputPath,
                               const Array& specificIntensityv, double n)
{
	const string tab = "\t";
	const int samples = 200;

	double T = 10;
	double factor = pow(1000000. / 10., 1. / samples);

	ofstream output = IOTools::ofstreamFile(outputPath + "heatingcurve.dat");
	output << "# 0temperature 1net 2heat 3cool"
	       << " 4lineNet 5lineHeat 6lineCool"
	       << " 7continuumNet 8continuumHeat 9continuumCool"
	       << " 10ionizedFrac" << endl;
	for (int N = 0; N < samples; N++, T *= factor)
	{
		GasInterfaceImpl::Solution s = gi.calculateDensities(n, T, specificIntensityv);
		double heat = gi.heating(s);
		double cool = gi.cooling(s);
		double lHeat = gi.lineHeating(s);
		double lCool = gi.lineCooling(s);
		double cHeat = gi.continuumHeating(s);
		double cCool = gi.continuumCooling(s);

		double netHeating = heat - cool;
		double netLine = lHeat - lCool;
		double netCont = cHeat - cCool;

		output << T << tab << netHeating << tab << heat << tab << cool << tab << netLine
		       << tab << lHeat << tab << lCool << tab << netCont << tab << cHeat << tab
		       << cCool << tab << s.f << endl;
	}
	output.close();

	double isrf = TemplatedUtils::integrate<double>(gi.frequencyv(), specificIntensityv);

	cout << "Calculated heating curve under isrf of " << isrf << " erg / s / cm2 / sr = "
	     << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl;
}
