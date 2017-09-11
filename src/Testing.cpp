#include "Testing.h"
#include "ChemicalNetwork.h"
#include "ChemistrySolver.h"
#include "Constants.h"
#include "EigenAliases.h"
#include "Error.h"
#include "FreeBound.h"
#include "GasInterface.h"
#include "GasInterfaceImpl.h"
#include "GrainInterface.h"
#include "GrainPhotoelectricEffect.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "NLevel.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"

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

Array Testing::freqToWavGrid(const Array& frequencyv)
{
	size_t numWav = frequencyv.size();
	Array wavelengthv(numWav);
	for (size_t iWav = 0; iWav < numWav; iWav++)
	{
		wavelengthv[iWav] = Constant::LIGHT / frequencyv[numWav - 1 - iWav];
	}
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

		/* Skip the wing points of line if it lies near the core of the previous one, to
		   prevent too much points from bunching up. */
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

Array Testing::generateSpecificIntensityv(const Array& frequencyv, double Tc, double G0)
{
	Array I_nu(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		I_nu[iFreq] = SpecialFunctions::planck(frequencyv[iFreq], Tc);

	// Cut out the UV part
	size_t i = 0;
	size_t startUV, endUV;
	while (frequencyv[i] < Constant::LIGHT / (2400 * Constant::ANG_CM) && i < frequencyv.size())
		i++;
	startUV = i > 0 ? i - 1 : 0;
	while (frequencyv[i] < Constant::LIGHT / (912 * Constant::ANG_CM) && i < frequencyv.size())
		i++;
	endUV = i + 1;
	cout << "UV goes from " << startUV << " to " << endUV << endl;
	vector<double> frequenciesUV(begin(frequencyv) + startUV, begin(frequencyv) + endUV);
	vector<double> isrfUV(begin(I_nu) + startUV, begin(I_nu) + endUV);

	// Integrate over the UV only
	double UVdensity = Constant::FPI / Constant::LIGHT *
	                   TemplatedUtils::integrate<double>(frequenciesUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING;

	// Rescale to _G0
	I_nu *= G0 / currentG0;

	vector<double> isrfUVbis(begin(I_nu) + startUV, begin(I_nu) + endUV);

	// Integrate over the UV only
	double UVdensitybis = Constant::FPI / Constant::LIGHT *
	                      TemplatedUtils::integrate<double>(frequenciesUV, isrfUVbis);
	cout << "Normalized spectrum uv = " << UVdensitybis << " ("
	     << UVdensitybis / Constant::HABING << " habing)" << endl;

	// Write out the ISRF
	ofstream out = IOTools::ofstreamFile("testing/isrf.txt");
	for (size_t b = 0; b < I_nu.size(); b++)
		out << frequencyv[b] << '\t' << I_nu[b] << '\n';
	out.close();
	out = IOTools::ofstreamFile("testing/isrfUV.txt");
	for (size_t b = 0; b < isrfUV.size(); b++)
		out << frequenciesUV[b] << '\t' << isrfUVbis[b] << '\n';
	out.close();
	return I_nu;
}

Array Testing::freqToWavSpecificIntensity(const Array& frequencyv,
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
	double n = 1;
	double f = 1;
	vector<double> kT_eVv = generateGeometricGridv(300, 1e-3, 1e3);
	for (double kT_eV : kT_eVv)
	{
		double T = kT_eV / Constant::BOLTZMAN / Constant::ERG_EV;
		double cool = Ionization::cooling(n * (1 - f), f * n, f * n, T) / Constant::RYDBERG;
		out << kT_eV << "\t" << cool << endl;
	}
	out.close();
}

void Testing::runGasInterfaceImpl(const GasModule::GasInterface& gi, const std::string& outputPath,
                                  double Tc, double G0, double n, double expectedTemperature)
{
	const Array& frequencyv = gi.frequencyv();

	Array specificIntensityv = generateSpecificIntensityv(frequencyv, Tc, G0);

	GasModule::GasState gs;
	GasModule::GrainInterface grainInfo{};
	gi.updateGasState(gs, n, expectedTemperature, specificIntensityv, grainInfo);

	cout << "Equilibrium temperature: " << gs.temperature() << endl;
	cout << "Ionized fraction: " << gs.ionizedFraction() << endl;

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

	cout << "----------------------------------" << endl;
	cout << "plotting heating curve..." << endl;
	plotHeatingCurve(*gi.pimpl(), outputPath, specificIntensityv, n);
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
		       << cCool << tab << gi.f(s) << endl;
	}
	output.close();

	double isrf = TemplatedUtils::integrate<double>(gi.frequencyv(), specificIntensityv);

	cout << "Calculated heating curve under isrf of " << isrf << " erg / s / cm2 / sr = "
	     << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl;
}

void Testing::testPhotoelectricHeating()
{
	double n = 2.5e1;
	double f = 3.e-4;
	double ne = n * f;
	double gasT = 1000;
	vector<double> G0values;
	if (gasT == 1000)
		G0values = {2.45e-2, 2.45e-1, 2.45e0, 2.45e1, 2.45e2};
	if (gasT == 100)
		G0values = {.75e-1, .75e0, .75e1, .75e2, .75e3};

	GrainPhotoelectricEffect phr{GasModule::GrainType::CAR};
	phr.yieldFunctionTest();
	for (double G0 : G0values)
	{
		phr.chargeBalanceTest(G0, gasT, ne, ne);
		phr.heatingRateTest(G0, gasT, ne);
	}
}

void Testing::testACollapse()
{
	HydrogenFromFiles hff(0);
	EMatrix avv = hff.avv().array();
	cout << hff.avv() << endl;
	cout << "Compare this with values directly from NIST below:" << endl;
	EMatrix nistA(5, 5);
	// clang-format off
	nistA << 0, 0, 0, 0, 0, 4.6986e+08, 0, 0, 0, 0, 5.5751e+07, 4.4101e+07, 0, 0, 0, 1.2785e+07, 8.4193e+06, 8.9860e+06, 0, 0, 4.1250e+06, 2.5304e+06, 2.2008e+06, 2.6993e+06, 0;
	// clang-format on
	cout << nistA << endl;
	cout << "The element-wise relative difference is " << endl;
	EMatrix relDiff = (avv - nistA).array() / nistA.array();
	// Take out all the nan's
	for (int i = 0; i < relDiff.size(); i++)
	{
		auto* pointer = relDiff.data() + i;
		auto value = *pointer;
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
	EMatrix avv = hff.avv();
	EMatrix cvv = hff.cvv(T, ne, np);

	/* Calculate and write out (q_n(l-1) + q_n(l+1)) / A_nl, where A_nl is the total downwards
	 rate from level nl. */
	EVector anlv = avv.rowwise().sum();
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

void Testing::testChemistry()
{
	const double T = 10000;
	vector<double> freqvec = generateGeometricGridv(200, 1e11, 1e16);
	Array frequencyv(freqvec.data(), freqvec.size());
	Array specificIntensityv = generateSpecificIntensityv(frequencyv, 25000, 10);

	ChemistrySolver cs(make_unique<ChemicalNetwork>());

	// Apply an artificial dissociation. If chemistry is correct, all H2 should be converted
	// into H.
	double kdiss = 0;
	EVector kv = cs.chemicalNetwork()->rateCoeffv(T, frequencyv, specificIntensityv, kdiss);
	cout << "Rate coeff: ionization, recombination, dissociation" << endl << kv << endl;

	int ie = ChemicalNetwork::speciesIndexm.at("e-");
	int ip = ChemicalNetwork::speciesIndexm.at("H+");
	int iH = ChemicalNetwork::speciesIndexm.at("H");
	int iH2 = ChemicalNetwork::speciesIndexm.at("H2");

	EVector n0v(4);
	n0v(ie) = 0;
	n0v(ip) = 0;
	n0v(iH) = 100;
	n0v(iH2) = 0;

	EVector nv = cs.solveBalance(kv, n0v);
	double ionizedFraction = Ionization::solveBalance(nv(iH) + nv(ip), T, frequencyv,
	                                                  specificIntensityv);

	cout << "Compare with ionized fraction calculation: " << endl;
	cout << "f = " << ionizedFraction << endl;
	assert(abs(nv(ip) / (nv(ip) + nv(iH)) - ionizedFraction) < 0.01);
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
	EVector evhc = hhc.ev();
	EVector evff = hff.ev();
	hc_vs_ff(evhc, evff);

	assert(hhc.gv() == hff.gv());

	cout << "A coefficients:" << endl;
	EMatrix avvhc = hhc.avv();
	EMatrix avvff = hff.avv();
	hc_vs_ff(avvhc, avvff);

	cout << "Extra A:" << endl;
	EMatrix eavvhc = hhc.extraAvv();
	EMatrix eavvff = hff.extraAvv();
	hc_vs_ff(eavvhc, eavvff);

	double T = 1e4;
	double ne = 1e4;
	double np = 1e4;

	cout << "Collisions:" << endl;
	EMatrix cvvhc = hhc.cvv(T, ne, np);
	EMatrix cvvff = hff.cvv(T, ne, np);
	hc_vs_ff(cvvhc, cvvff);

	cout << "Recombinations:" << endl;
	EVector alphavhc = hhc.sourcev(T, ne, np) / ne / np;
	EVector alphavff = hff.sourcev(T, ne, np) / ne / np;
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

	GasModule::GasInterface gihhc(frequencyv, "hhc", false);
	runGasInterfaceImpl(gihhc, "hardcoded/");

	GasModule::GasInterface gihff(frequencyv, "hff2", false);
	runGasInterfaceImpl(gihff, "fromfiles/");
}

GasModule::GasInterface Testing::genFullModel()
{
	bool molecular = true;

	vector<double> tempFrequencyv =
	                generateGeometricGridv(2000, Constant::LIGHT / (1e4 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.005 * Constant::UM_CM));
	Array unrefined(tempFrequencyv.data(), tempFrequencyv.size());

	HydrogenLevels hl(make_shared<HydrogenFromFiles>(), unrefined);
	FreeBound fb(unrefined);
	Array frequencyv = improveFrequencyGrid(hl, fb, unrefined);

	return {frequencyv, "", molecular};
}

void Testing::runFullModel() { runGasInterfaceImpl(genFullModel(), ""); }

void Testing::runWithDust()
{
	GasModule::GasInterface gasInterface{genFullModel()};
	double nHtotal = 10;

	std::vector<GasModule::GrainInterface::Population> grainPopv;

	Array sizev, densityv, temperaturev;
	// Provide sizes in cm
	sizev = {1e-7, 1e-6, 1e-5};
	densityv = {0.01, 0.005, 0.0025}; // number of grain per H atom
	densityv *= nHtotal; // density in cm-3
	// Actually, a temperature distribution per grain size is the most detailed form we'll be
	// using. Maybe we need multiple versions of the grain interface with regards to storing the
	// grain temperatures.
	temperaturev = {50, 20, 10};
	// And now I need some absorption efficiencies for every wavelength. Let's try to use the
	// old photoelectric heating test code.
	std::vector<Array> qAbsvv{GrainPhotoelectricEffect::qAbsvvForTesting(
	                sizev, gasInterface.frequencyv())};

	grainPopv.emplace_back(GasModule::GrainType::CAR, sizev, densityv, temperaturev, qAbsvv);

	GasModule::GrainInterface grainInterface{grainPopv};

	// TODO: generate a radiation field and run that thing
}
