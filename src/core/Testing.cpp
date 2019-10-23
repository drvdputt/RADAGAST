#include "Testing.hpp"
#include "BigH2Model.hpp"
#include "FreeBound.hpp"
#include "GasDiagnostics.hpp"
#include "GasInterface.hpp"
#include "GasInterfaceImpl.hpp"
#include "GasStruct.hpp"
#include "GrainPhotoelectricEffect.hpp"
#include "GrainType.hpp"
#include "HFromFiles.hpp"
#include "IOTools.hpp"
#include "Ionization.hpp"
#include "SpecialFunctions.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"

#include <gsl/gsl_const_cgs.h>

#include <sys/stat.h>

using namespace std;

namespace
{
typedef struct QabsDataSet
{
	vector<double> FILELAMBDAV;
	vector<double> FILEAV;
	vector<vector<double>> QABSVV;

} QabsDataSet;

QabsDataSet readQabs(bool car)
{
	QabsDataSet data;

	// We dont need these, but store them anyway for clarity
	vector<vector<double>> QSCAVV, ASYMMPARVV;

	///////////////////// Begin copy-paste from SKIRT
	bool reverse = true;
	bool skip1 = false;
	bool skip2 = false;
	bool skip3 = false;
	size_t _Nlambda, _Na;
	// open the file
	ifstream file;
	if (car)
		file = IOTools::ifstreamRepoFile("dat/Gra_81.dat");
	else
		file = IOTools::ifstreamRepoFile("dat/suvSil_81.dat");

	// skip header lines and read the grid size
	string line;
	while (file.peek() == '#')
		getline(file, line);
	file >> _Na;
	getline(file, line); // ignore anything else on this line
	file >> _Nlambda;
	getline(file, line); // ignore anything else on this line

	// resize the vectors
	data.FILELAMBDAV.resize(_Nlambda);
	data.FILEAV.resize(_Na);
	data.QABSVV.resize(_Nlambda, vector<double>(_Na));
	QSCAVV.resize(_Nlambda, vector<double>(_Na));
	ASYMMPARVV.resize(_Nlambda, vector<double>(_Na));

	// determine the loop conditions for wavelength lines
	int kbeg = reverse ? _Nlambda - 1 : 0;
	int kend = reverse ? -1 : _Nlambda;
	int kinc = reverse ? -1 : 1;

	// read the data blocks
	double dummy;
	for (size_t i = 0; i < _Na; i++)
	{
		file >> data.FILEAV[i];
		data.FILEAV[i] *= 1e-6; // convert from micron to m
		getline(file, line); // ignore anything else on this line

		for (int k = kbeg; k != kend; k += kinc)
		{
			if (skip1)
				file >> dummy;
			file >> data.FILELAMBDAV[k];
			data.FILELAMBDAV[k] *= 1e-6; // convert from micron to m
			if (skip2)
				file >> dummy;
			file >> data.QABSVV[k][i];
			file >> QSCAVV[k][i];
			if (skip3)
				file >> dummy;
			file >> ASYMMPARVV[k][i];
			getline(file, line); // ignore anything else on this line
		}
	}

	// close the file
	file.close();
	///////////////////// End copy-paste from SKIRT

	// Convert the wavelengths and grain sizes from microns to centimeters
	for (double& d : data.FILELAMBDAV)
		d *= 100.; // m to cm
	for (double& d : data.FILEAV)
		d *= 100.; // m to cm

	return data;
}

Array generateQabsv(const QabsDataSet& qds, double a, const Array& frequencyv)
{
	// Annoying edge case included
	if (a <= qds.FILEAV[0])
		return Array(frequencyv.size());

	if (a > qds.FILEAV[qds.FILEAV.size() - 1])
		return Array(frequencyv.size());

	size_t aIndex = TemplatedUtils::index(a, qds.FILEAV);
	assert(aIndex > 0);

	Array wavelengthv = Testing::freqToWavGrid(frequencyv);
	vector<double> QabsWav(wavelengthv.size());
	for (size_t w = 0; w < wavelengthv.size(); w++)
	{
		double wav = wavelengthv[w];

		size_t wIndex = TemplatedUtils::index(wav, qds.FILELAMBDAV);
		double wLeft = qds.FILELAMBDAV[wIndex - 1];
		double wRight = qds.FILELAMBDAV[wIndex];

		double aLow = qds.FILEAV[aIndex - 1];
		double aUp = qds.FILEAV[aIndex];

		const auto& q = qds.QABSVV;
		QabsWav[w] = TemplatedUtils::interpolateRectangular(
		                wav, a, wLeft, wRight, aLow, aUp, q[wIndex - 1][aIndex - 1],
		                q[wIndex][aIndex - 1], q[wIndex - 1][aIndex],
		                q[wIndex][aIndex]);
	}

// #define PLOT_QABS
#ifdef PLOT_QABS
	stringstream filename;
	filename << "photoelectric/multi-qabs/qabs_a" << setfill('0') << setw(8)
	         << setprecision(2) << fixed << a / Constant::ANG_CM << ".txt";
	ofstream qabsfile = IOTools::ofstreamFile(filename.str());
	for (size_t i = 0; i < frequencyv.size(); i++)
		qabsfile << frequencyv[i] * Constant::CM_UM << '\t' << QabsWav[i] << endl;
	qabsfile.close();
#endif
	// Reverse order to make Qabs a function of frequency index
	reverse(begin(QabsWav), end(QabsWav));
	return Array(QabsWav.data(), QabsWav.size());
}

} // namespace

vector<Array> Testing::qAbsvvForTesting(bool car, const Array& av, const Array& frequencyv)
{
	// Choose carbon
	QabsDataSet qds = readQabs(car);

	vector<Array> result;
	result.reserve(av.size());
	for (double a : av)
		result.emplace_back(generateQabsv(qds, a, frequencyv));
	return result;
}

double Testing::equilibriumTemperature(const Array& frequencyv, const Array& specificIntensityv,
                                       const Array& crossSectionv)
{
	// Assumes that both arguments on the same frequency grid

	// There's probably quite some mistakes in this (instead of finding the right way to
	// integrate / averaging, I just winged it). But hopefully it's close enough to make the
	// cloudy comparison easier.

	double absorption = TemplatedUtils::integrate<double, Array, Array>(
	                frequencyv, crossSectionv * specificIntensityv);

	auto heating = [&](double T) -> int {
		Array blackbodyIntegrandv(frequencyv.size());
		for (size_t i = 0; i < frequencyv.size(); i++)
			blackbodyIntegrandv[i] = SpecialFunctions::planck(frequencyv[i], T) *
			                         crossSectionv[i];

		double bbEmission = 0;
		if (T > 0.)
			bbEmission = TemplatedUtils::integrate<double, Array, Array>(
			                frequencyv, blackbodyIntegrandv);
		if (bbEmission < absorption)
			return 1;
		if (bbEmission > absorption)
			return -1;
		else
			return 0;
	};

	double result = TemplatedUtils::binaryIntervalSearch<double>(heating, 30., 1.e-3, 300.,
	                                                             1.);
	return result;
}

Array Testing::improveFrequencyGrid(const LevelCoefficients& boundBound, const Array& oldPoints)
{
	// Add extra points for the lines
	int numLines;
	Array lineFreqv, naturalWidthv;
	boundBound.lineInfo(numLines, lineFreqv, naturalWidthv);

	double lineWindowFactor = 1.;
	double thermalFactor = sqrt(Constant::BOLTZMAN * 50000 / Constant::HMASS_CGS) /
	                       Constant::LIGHT;
	naturalWidthv = lineWindowFactor * (naturalWidthv + lineFreqv * thermalFactor);

	vector<double> gridVector(begin(oldPoints), end(oldPoints));
	Testing::refineFrequencyGrid(gridVector, 13, 2.5, lineFreqv, naturalWidthv);

	return Array(gridVector.data(), gridVector.size());
}

Array Testing::improveFrequencyGrid(const FreeBound& freeBound, const Array& oldPoints)
{
	// Copy the current points
	vector<double> gridVector(begin(oldPoints), end(oldPoints));

	// Add points for the jumps in the bound-bound spectrum
	const Array& thresholdv = freeBound.thresholdv();
	// Don't bother with the last jump, because that's the end of the data
	Array jumpFreqv(&thresholdv[0], thresholdv.size() - 2);
	Testing::refineFrequencyGrid(gridVector, 3, 1., jumpFreqv, 1e-6 * jumpFreqv);

	// And for the ionization threshold
	Array ionThr({Ionization::THRESHOLD});
	Testing::refineFrequencyGrid(gridVector, 3, 1., ionThr, 1e-6 * ionThr);

	return Array(gridVector.data(), gridVector.size());
}

Array Testing::generateGeometricGridv(size_t nPoints, double min, double max)
{
	Array frequencyv(nPoints);
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

Array Testing::defaultFrequencyv(size_t numPoints)
{
	Array coarsev = generateGeometricGridv(numPoints, defaultMinFreq, defaultMaxFreq);
	return coarsev;
}

Array Testing::freqToWavGrid(const Array& frequencyv)
{
	size_t numWav = frequencyv.size();
	Array wavelengthv(numWav);
	for (size_t iWav = 0; iWav < numWav; iWav++)
		wavelengthv[iWav] = Constant::LIGHT / frequencyv[numWav - 1 - iWav];
	return wavelengthv;
}

void Testing::refineFrequencyGrid(vector<double>& grid, size_t nPerLine, double spacingPower,
                                  Array lineFreqv, Array freqWidthv)
{
	// Function that inserts a point, but only if the neighbouring points are not
	// too close (distance smaller than the resolution argument)
	auto insertIfRoom = [&](double freq, double resolution) {
		// Find the points that would neighbour our candidate frequency in the
		// current grid
		size_t iRight = TemplatedUtils::index(freq, grid);
		double fRight = grid[iRight];

		// If the point to the right is too close, do not insert
		if (abs(freq - fRight) < resolution)
			return;

		// Consider the point to the left, if there is one
		if (iRight > 0)
		{
			size_t iLeft = iRight - 1;
			double fLeft = grid[iLeft];
			if (abs(freq - fLeft) < resolution)
				return;
		}

		// We have now covered all exceptions (left too close, right too close)
		// and decide to insert the point.
		TemplatedUtils::sortedInsert(freq, grid);
	};

	// We want an odd nPerLine, so we can put 1 point in the center
	if (!(nPerLine % 2))
		nPerLine += 1;

	for (size_t i = 0; i < lineFreqv.size(); i++)
	{
		// Add a point at the center of the line, while keeping the vector sorted
		TemplatedUtils::sortedInsert<double, vector<double>>(lineFreqv[i], grid);

		if (nPerLine > 1)
		{
			double center = lineFreqv[i];
			double width = freqWidthv[i];
			// do not add points if there are already points at this distance
			double resolutionLimit = 0.1 * width / nPerLine;

			// Add the rest of the points in a power law spaced way
			size_t nOneSide = (nPerLine - 1) / 2;

			// Normalization factor. Last point is placed at distance 'width'
			double a = width / pow(nOneSide, spacingPower);
			for (size_t sidePoint = 1; sidePoint <= nOneSide; sidePoint++)
			{
				double distance = a * pow(sidePoint, spacingPower);
				// Left and right of center
				insertIfRoom(center - distance, resolutionLimit);
				insertIfRoom(center + distance, resolutionLimit);
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
	while (frequencyv[i] < Constant::LIGHT / (2400 * Constant::ANG_CM) &&
	       i < frequencyv.size())
		i++;
	startUV = i > 0 ? i - 1 : 0;
	while (frequencyv[i] < Constant::LIGHT / (912 * Constant::ANG_CM) &&
	       i < frequencyv.size())
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
	// ofstream out = IOTools::ofstreamFile("testing/isrf.txt");
	// for (size_t b = 0; b < I_nu.size(); b++)
	//	out << frequencyv[b] << '\t' << I_nu[b] << '\n';
	// out.close();
	// out = IOTools::ofstreamFile("testing/isrfUV.txt");
	// for (size_t b = 0; b < isrfUV.size(); b++)
	//	out << frequenciesUV[b] << '\t' << isrfUVbis[b] << '\n';
	// out.close();
	return I_nu;
}

Array Testing::freqToWavSpecificIntensity(const Array& frequencyv,
                                          const Array& specificIntensity_nu)
{
	Array I_lambda(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		I_lambda[I_lambda.size() - iFreq - 1] = specificIntensity_nu[iFreq] *
		                                        frequencyv[iFreq] * frequencyv[iFreq] /
		                                        Constant::LIGHT;
	return I_lambda;
}

void Testing::plotIonizationStuff()
{
	ofstream out;

	double A0 = 6.30e-18;
	out = IOTools::ofstreamFile("ionization/crossSection.dat");
	double tr = Ionization::THRESHOLD;
	for (double freq = 0; freq < 6.e16; freq += 1e14)
	{
		double sigma = Ionization::crossSection(freq);
		double eps = sqrt(freq / tr - 1);
		double sigmaTheoretical;
		if (freq > tr)
			sigmaTheoretical = A0 * pow(tr / freq, 4.) *
			                   exp(4 - 4 * atan(eps) / eps) /
			                   (1 - exp(-2 * Constant::PI / eps));
		else
			sigmaTheoretical = 0;
		out << freq << "\t" << sigma << "\t" << sigmaTheoretical << "\t"
		    << (sigma - sigmaTheoretical) / sigma << endl;
	}
	out.close();

	out = IOTools::ofstreamFile("ionization/recombinationCooling.dat");
	out << "#kT / eV \t alpha" << endl;
	// by choosing these parameters, panel 2 of figure 9 from 2017-Mao should be reproduced
	double n = 1;
	double f = 1;
	Array kT_eVv = generateGeometricGridv(300, 1e-3, 1e3);
	for (double kT_eV : kT_eVv)
	{
		double T = kT_eV / Constant::BOLTZMAN / Constant::ERG_EV;
		double cool = Ionization::cooling(n * (1 - f), f * n, f * n, T) /
		              Constant::RYDBERG;
		out << kT_eV << "\t" << cool << endl;
	}
	out.close();
}

void Testing::runGasInterfaceImpl(const GasModule::GasInterface& gi,
                                  const std::string& outputPath, double Tc, double G0, double n)
{
	Array specificIntensityv = generateSpecificIntensityv(gi.iFrequencyv(), Tc, G0);

	GasModule::GasState gs;
	GasModule::GrainInterface gri{};
	gi.updateGasState(gs, n, specificIntensityv, gri);
	writeGasState(outputPath, gi, gs);
}

void Testing::writeGasState(const string& outputPath, const GasModule::GasInterface& gi,
                            const GasModule::GasState& gs)
{
	cout << "Equilibrium temperature: " << gs.temperature() << endl;

	const Array& eFrequencyv = gi.eFrequencyv();
	const Array& emv = gs._emissivityv;
	Spectrum opv(gi.oFrequencyv(), gs._opacityv);

	char tab = '\t';
	ofstream out = IOTools::ofstreamFile(outputPath + "opticalProperties.dat");
	vector<std::string> colnames = {
	                "frequency",
	                "wavelength",
	                "intensity j_nu (erg s-1 cm-3 Hz-1 sr-1)",
	                "opacity alpha_nu (cm-1)",
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
		double freq = eFrequencyv[iFreq];
		double wav = Constant::LIGHT / freq * Constant::CM_UM;
		out.precision(9);
		out << scientific << freq << tab << wav << tab
		    << Constant::FPI * freq * emv[iFreq] << tab << opv.evaluate(freq) << endl;
		wavfile.precision(9);
		wavfile << wav << tab << freq << endl;
	}
	out.close();
	wavfile.close();

	out = IOTools::ofstreamFile(outputPath + "raw_opacity.dat");
	for (size_t iFreq = 0; iFreq < gs._opacityv.size(); iFreq++)
	{
		out << gi.oFrequencyv()[iFreq] << tab << gs._opacityv[iFreq] << endl;
	}
	out.close();

	// Print some line intensities relative to Hbeta
	double fHalpha = Constant::LIGHT / 656.453e-7;
	double fHbeta = Constant::LIGHT / 486.264e-7;
	double fHgamma = Constant::LIGHT / 434.165e-7;

	double fLya = Constant::LIGHT / 121.567e-7;

	double fPalpha = Constant::LIGHT / 1875.61e-7;
	double fPbeta = Constant::LIGHT / 1282.16e-7;

	double fBralpha = Constant::LIGHT / 4052.27e-7;

	auto evaluateSpectrum = [&](double f) {
		return TemplatedUtils::evaluateLinInterpf<double>(f, eFrequencyv, emv);
	};

	double Hbeta = evaluateSpectrum(fHbeta);
	cout << "Halpha / Hbeta " << evaluateSpectrum(fHalpha) / Hbeta << endl;
	cout << "Hgamma / Hbeta " << evaluateSpectrum(fHgamma) / Hbeta << endl;

	cout << "Lyalpha / Hbeta " << evaluateSpectrum(fLya) / Hbeta << endl;

	cout << "Palpha / Hbeta " << evaluateSpectrum(fPalpha) / Hbeta << endl;
	cout << "Pbeta / Hbeta " << evaluateSpectrum(fPbeta) / Hbeta << endl;

	cout << "Bralpha / HBeta " << evaluateSpectrum(fBralpha) / Hbeta << endl;
}

void Testing::writeGrains(const std::string& outputPath, const GasModule::GrainInterface& gr,
                          bool bulkCar)
{
	// Assume a certain bulk density (values taken from SKIRT source code for Draine
	// Graphite and Draine Silicate)
	double bulkDen;
	if (bulkCar)
		bulkDen = 2.24; // g cm-3
	else
		bulkDen = 3.0;

	// write out grain size distribution, preferable in g cm-3
	for (size_t i = 0; i < gr.numPopulations(); i++)
	{
		ColumnFile f(outputPath + "grainpop_" + std::to_string(i) + ".dat",
		             {"size", "density(cm-3)", "massdensity(g cm-3)", "temperature"});

		auto pop = gr.population(i);
		for (size_t m = 0; m < pop->numSizes(); m++)
		{
			double a = pop->size(m);
			double numberDen = pop->density(m);
			double mass = 4. / 3. * Constant::PI * a * a * a * bulkDen;
			f.writeLine<Array>(
			                {a, numberDen, numberDen * mass, pop->temperature(m)});
		}
	}
}

void Testing::plotHeatingCurve_main()
{
	Array frequencyv =
	                generateGeometricGridv(500, Constant::LIGHT / (1e3 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.005 * Constant::UM_CM));
	double Tc = 30000;
	double g0 = 1e0;
	double n = 1000;
	Spectrum specificIntensity{frequencyv, generateSpecificIntensityv(frequencyv, Tc, g0)};
	GasModule::GasInterface gi{frequencyv, frequencyv, frequencyv, "", "none"};
	GasModule::GrainInterface gri{};
	string outputPath = "heatingcurve/";
	plotHeatingCurve(*gi.pimpl(), outputPath, n, specificIntensity, gri);
}

void Testing::plotHeatingCurve(const GasInterfaceImpl& gi, const std::string& outputPath,
                               double n, const Spectrum& specificIntensity,
                               GasModule::GrainInterface& gri)
{
	// Initial
	Array Tv = Testing::generateGeometricGridv(100, 10, 50000);
	double T0 = 10000;
	GasSolution s = gi.solveDensities(n, T0, specificIntensity, gri);
	GasDiagnostics gd;
	s.fillDiagnostics(&gd);

	ColumnFile heatFile(outputPath + "heatcool.dat",
	                    {"temperature", "net", "heat", "cool", "grainheat"});
	ColumnFile densFile(outputPath + "densities.dat",
	                    {"temperature", "e-", "H+", "H", "H2", "Htot"});

	vector<string> rateFileColNames = {"temperature"};
	for (auto name : gd.reactionNames())
		rateFileColNames.emplace_back(name);
	ColumnFile rateFile(outputPath + "chemrates.dat", rateFileColNames);

	std::vector<int> h_conservation = {0, 1, 1, 2};
	size_t numSpecies = h_conservation.size();
	vector<double> densFileLine(2 + numSpecies);

	const GasSolution* previous = nullptr;
	auto outputCooling = [&](double t) {
		std::cout << "T = " << t << '\n';
		gi.solveDensities(s, n, t, specificIntensity, gri, previous);
		s.fillDiagnostics(&gd);

		double heat = s.heating();
		double cool = s.cooling();
		double grainHeat = gd.photoelectricHeating().sum();
		double netHeating = heat - cool + grainHeat;
		heatFile.writeLine<Array>({t, netHeating, heat, cool, grainHeat});

		densFileLine[0] = t;
		double totalH = 0;
		for (int i = 0; i < numSpecies; i++)
		{
			densFileLine[1 + i] = s.speciesNv()(i);
			totalH += h_conservation[i] * s.speciesNv()(i);
		}
		densFileLine[1 + numSpecies] = totalH;
		densFile.writeLine(densFileLine);

		vector<double> lineValues;
		lineValues.reserve(rateFileColNames.size());
		lineValues.emplace_back(t);
		for (double rate : gd.reactionRates())
			lineValues.emplace_back(rate);
		rateFile.writeLine(lineValues);
	};

	// forward sweep
	for (double& d : Tv)
		outputCooling(d);
	// reverse sweep
	size_t i = Tv.size();
	while (i--)
		outputCooling(Tv[i]);
}

void Testing::plotPhotoelectricHeating()
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
	vector<int> pickValues{0, 1, 2, 3, 4};

	unique_ptr<GrainType> grainType{
	                GrainTypeFactory::makeBuiltin(GasModule::GrainTypeLabel::CAR)};
	GrainPhotoelectricEffect phr(*grainType);
	phr.yieldFunctionTest();

	for (int i : pickValues)
	{
		phr.chargeBalanceTest(G0values[i], gasT, ne, ne);
		phr.heatingRateTest(G0values[i], gasT, ne);
	}
}

void Testing::plotInterpolationTests()
{
	double mn = 1;
	double mx = 10;

	Array baseGridv = generateGeometricGridv(30, mn, mx);
	Array baseFv(baseGridv.size());
	for (size_t i = 0; i < baseGridv.size(); i++)
		baseFv[i] = exp(-baseGridv[i] * baseGridv[i]);

	auto out = IOTools::ofstreamFile("base.dat");
	for (size_t i = 0; i < baseGridv.size(); i++)
		out << baseGridv[i] << '\t' << baseFv[i] << endl;
	out.close();

	Spectrum s(baseGridv, baseFv);

	Array finerGridv = generateGeometricGridv(1000, mn, mx);
	Array finerFv = s.binned(finerGridv);

	out = IOTools::ofstreamFile("finer.dat");
	for (size_t i = 0; i < finerGridv.size(); i++)
		// out << finerGridv[i] << '\t' << finerFv[i] << endl;
		out << finerGridv[i] << '\t' << s.evaluate(finerGridv[i]) << endl;
	out.close();

	Array coarserGridv = generateGeometricGridv(10, mn, mx);
	Array coarserFv = s.binned(coarserGridv);

	out = IOTools::ofstreamFile("coarser.dat");
	for (size_t i = 0; i < coarserGridv.size(); i++)
		out << coarserGridv[i] << '\t' << coarserFv[i] << endl;
	out.close();
}

void Testing::plotPS64Collisions()
{
	const double T = 10000;
	const double ne = 1e4;
	const double np = 1e4;

	HFromFiles hff(5);
	EMatrix avv = hff.avv();
	EVector speciesNv = EVector::Zero(SpeciesIndex::size());
	speciesNv(SpeciesIndex::ine()) = ne;
	speciesNv(SpeciesIndex::inp()) = np;
	GasStruct gas(T, speciesNv);
	EMatrix cvv = hff.cvv(gas);

	/* Calculate and write out (q_n(l-1) + q_n(l+1)) / A_nl, where A_nl is the total
	   downwards rate from level nl. */
	EVector anlv = avv.rowwise().sum();
	for (int n : std::array<int, 2>{4, 5})
	{
		ofstream out = IOTools::ofstreamFile("ps64/t" + to_string(n) + "l_q" +
		                                     to_string(n) + "l.dat");
		for (int li = 0; li < n; li++)
		{
			size_t nliIndex = hff.index(n, li);
			// Decay rate to all other levels
			double anl = anlv(nliIndex);
			double qnl = 0;

			// Get the total decay rate due to l-changing collisions
			for (int lf = 0; lf < n; lf++)
			{
				size_t nlfIndex = hff.index(n, lf);
				qnl += cvv(nliIndex, nlfIndex);

				// l-changing collision rates should be zero except for changes
				// by 1
				int deltal = li - lf;
				if (cvv(nliIndex, nlfIndex) > 0 && abs(deltal) != 1)
				{
					stringstream ss;
					ss << "The l-changing coefficient from " << n << ","
					   << li << " to " << n << "," << lf
					   << " should be zero.";
					Error::runtime(ss.str());
				}
			}
			out << li << "\t" << qnl / anl << "\t" << qnl << "\t" << anl << endl;
		}
		out.close();
	}
}

void Testing::runH2(bool write)
{
	cout << "RUN_H2" << endl;

	int maxJ = 10;
	int maxV = 3;

	double nH2 = 1000;
	double ne = 100;
	double np = 100;
	double nH = 100;
	double T = 500;
	double Tc = 10000;
	double G0 = 10;

	// Base grid
	Array unrefinedv =
	                generateGeometricGridv(20000, Constant::LIGHT / (1e4 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.005 * Constant::UM_CM));

	Array specificIntensityv = generateSpecificIntensityv(unrefinedv, Tc, G0);
	Spectrum specificIntensity(unrefinedv, specificIntensityv);

	// Add points for H2 lines
	H2Data h2l(maxJ, maxV);
	Array frequencyv = improveFrequencyGrid(h2l, unrefinedv);

	// Set the densities
	EVector speciesNv{EVector::Zero(SpeciesIndex::size())};
	speciesNv(SpeciesIndex::inH2()) = nH2;
	speciesNv(SpeciesIndex::ine()) = ne;
	speciesNv(SpeciesIndex::inp()) = np;
	speciesNv(SpeciesIndex::inH()) = nH;

	GasStruct gas(T, speciesNv);
	BigH2Model h2m(&h2l);
	h2m.solve(nH2, gas, specificIntensity);
	if (write)
	{
		Array emissivityv = h2m.emissivityv(frequencyv);
		Array opacityv = h2m.opacityv(unrefinedv);
		ofstream h2optical = IOTools::ofstreamFile("h2/opticalProperties.dat");
		h2optical << "# nu (hz)\tlambda (micron)\temissivity\topacity\tlineOpacity\n";
		const string tab{"\t"};
		for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		{
			h2optical << frequencyv[iFreq] << tab
			          << Constant::LIGHT / frequencyv[iFreq] * Constant::CM_UM
			          << tab << emissivityv[iFreq] << tab << opacityv[iFreq]
			          << '\n';
		}
	}
}

void Testing::runFromFilesvsHardCoded()
{
	Array unrefinedv =
	                generateGeometricGridv(1000, Constant::LIGHT / (1e10 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.00001 * Constant::UM_CM));

	HFromFiles hl(5);
	FreeBound fb;
	Array frequencyv = improveFrequencyGrid(hl, unrefinedv);
	frequencyv = improveFrequencyGrid(fb, frequencyv);

	GasModule::GasInterface gihhc(unrefinedv, unrefinedv, frequencyv, "hhc", "none");
	runGasInterfaceImpl(gihhc, "hardcoded/");

	GasModule::GasInterface gihff(unrefinedv, unrefinedv, frequencyv, "hff2", "none");
	runGasInterfaceImpl(gihff, "fromfiles/");
}

GasModule::GasInterface Testing::genFullModel(bool refine)
{
	Array coarsev = defaultFrequencyv(3000);
	Array eFrequencyv = coarsev;

	if (refine)
	{
		cout << "Constructing model to help with refining frequency grid" << endl;
		HFromFiles hl;
		FreeBound fb;
		H2Data h2l(99, 99);
		Array frequencyv = improveFrequencyGrid(hl, coarsev);
		frequencyv = improveFrequencyGrid(fb, frequencyv);
		frequencyv = improveFrequencyGrid(h2l, frequencyv);
	}

	cout << "Constructing new model using the improved frequency grid" << endl;
	return {coarsev, coarsev, eFrequencyv, "", "99 99"};
}

GasModule::GasInterface Testing::genHonlyModel()
{
	Array coarsev = defaultFrequencyv();

	cout << "Constructing H model to help with refining frequency grid" << endl;
	HFromFiles hl;
	FreeBound fb;

	Array frequencyv = improveFrequencyGrid(hl, coarsev);
	frequencyv = improveFrequencyGrid(fb, frequencyv);

	cout << "Constructing new model using the improved frequency grid" << endl;

	// TODO: make sure that this construction happens correctly. Since the rewrite, an Honly
	// model needs to be done differently.
	return {coarsev, coarsev, frequencyv, "", "none"};
}

void Testing::runFullModel()
{
	cout << "RUN_FULL_MODEL\n" << endl;
	GasModule::GasInterface gi = genHonlyModel();
	double Tc = 40000;
	double G0 = 100;
	double n = 1000;
	runGasInterfaceImpl(gi, "gasOnly/", Tc, G0, n);
}

GasModule::GrainInterface Testing::genMRNDust(double nHtotal, const Spectrum& specificIntensity)
{
	auto grainPopv{make_unique<vector<GasModule::GrainInterface::Population>>()};

	// need grains from .005 to .25 micron
	double amin = 50 * Constant::ANG_CM;
	double amax = 0.25 * Constant::UM_CM;

	// Shared properties
	size_t numSizes = 10;
	Array bin_edges = generateGeometricGridv(numSizes + 1, amin, amax);

	// Properties that differ between car and sil
	auto addMRNCarOrSil = [&](bool car) {
		double log10norm = car ? -25.13 : -25.11;
		double C = pow(10., log10norm);
		Array densityv(numSizes);
		Array sizev(numSizes);
		Array areav(numSizes);
		for (size_t i = 0; i < numSizes; i++)
		{
			// analytic integration of dn = C nH a^{-3.5} da
			// double bin_width = bin_edges[i + 1] - bin_edges[i];
			// densityv[i] = C * nHtotal * pow(sizev[i], -3.5) * bin_width;
			densityv[i] = C * nHtotal / 2.5 *
			              (pow(bin_edges[i], -2.5) - pow(bin_edges[i + 1], -2.5));

			auto average_power_of_size = [&](double power) {
				double exponent = 2.5 - power;
				return C * nHtotal / exponent *
				       (pow(bin_edges[i], -exponent) -
				        pow(bin_edges[i + 1], -exponent)) /
				       densityv[i];
			};
			sizev[i] = average_power_of_size(1.);
			areav[i] = average_power_of_size(2.);
		}
		auto label = car ? GasModule::GrainTypeLabel::CAR
		                 : GasModule::GrainTypeLabel::SIL;
		const auto& qabsvv =
		                qAbsvvForTesting(car, sizev, specificIntensity.frequencyv());

		Array temperaturev(numSizes);
		for (size_t i = 0; i < numSizes; i++)
		{
			Array crossSectionv = qabsvv[i] * Constant::PI * areav[i];
			temperaturev[i] = equilibriumTemperature(specificIntensity.frequencyv(),
			                                         specificIntensity.valuev(),
			                                         crossSectionv);
			cout << "grain " << i << ": " << temperaturev[i] << " K\n";
		}

		grainPopv->emplace_back(label, sizev, densityv, temperaturev, qabsvv);
	};

	// Carbonaceous population:
	addMRNCarOrSil(true);

	// Silicate population:
	// addMRNCarOrSil(false);

	GasModule::GrainInterface grainInterface(move(grainPopv));
	return grainInterface;
}

namespace
{
/** writes out a map as a single line column file, where the keys are the column names */
void writeMapAsColumnFile(string filename, const map<string, double>& map)
{
	vector<string> namev;
	vector<double> valuev;
	namev.reserve(map.size());
	valuev.reserve(map.size());
	for (auto const& item : map)
	{
		namev.emplace_back(item.first);
		valuev.emplace_back(item.second);
	}
	ColumnFile f(filename, namev);
	f.writeLine(valuev);
}
} // namespace

void Testing::runMRNDust(bool write, double nH, double Tc, double lumSol, bool own_dir)
{
	cout << "RUN_MRN_DUST\n";
	cout << "Parameters: nH " << nH << " Tc " << Tc << " lum " << lumSol << '\n';

	// Gas model
	GasModule::GasInterface gasInterface = genFullModel();
	double nHtotal = nH;

	// Radiation field
	// Tc argument is color temperature
	double bollum = lumSol * Constant::SOL_LUM;
	double distance = 1. * Constant::PC_CM; // 1.0 AU

	// double G0{1e2};
	Array frequencyv = gasInterface.iFrequencyv();
	Array specificIntensityv(frequencyv.size());
	for (int i = 0; i < frequencyv.size(); i++)
		specificIntensityv[i] = SpecialFunctions::planck(frequencyv[i], Tc);

	specificIntensityv *= bollum / 16 / Constant::PI / distance / distance /
	                      GSL_CONST_CGS_STEFAN_BOLTZMANN_CONSTANT / pow(Tc, 4);

	// Dust model
	auto gri = genMRNDust(nHtotal, Spectrum{frequencyv, specificIntensityv});

	auto gi_pimpl = gasInterface.pimpl();
	Spectrum specificIntensity(frequencyv, specificIntensityv);
	GasSolution s = gi_pimpl->solveTemperature(nHtotal, specificIntensity, gri);
	// Fixed temperature call, for convenience when testing:
	// GasSolution s = gi_pimpl->solveDensities(nHtotal, 49.4, specificIntensity, gri);
	if (write)
	{
		GasDiagnostics gd;
		s.fillDiagnostics(&gd);

		cout << "Htot = " << s.heating() << '\n';
		cout << "grainHeat = " << gd.photoelectricHeating().sum() << '\n';
		cout << "Ctot = " << s.cooling() << '\n';
		cout << "eden = " << s.ne() << '\n';
		cout << "H+ " << s.np() << '\n';
		cout << "HI " << s.nH() << '\n';
		cout << "H2 " << s.nH2() << '\n';
		for (size_t i = 0; i < gd.reactionNames().size(); i++)
			cout << gd.reactionNames()[i] << " = " << gd.reactionRates()[i] << '\n';

		string prefix = "MRNDust/";
		if (own_dir)
		{
			stringstream ss;
			ss << "nH" << nHtotal << "_Tc" << Tc << "_lum" << lumSol << '_'
			   << prefix;
			prefix = ss.str();
			if (prefix[prefix.size() - 1] == '/')
				mkdir(prefix.c_str(), 0755);
		}

		ColumnFile overview(prefix + "overview.dat", {"T", "eden", "H+", "HI", "H2"});
		overview.writeLine<Array>({s.t(), s.ne(), s.np(), s.nH(), s.nH2()});

		ColumnFile rates(prefix + "rates.dat", gd.reactionNames());
		rates.writeLine(gd.reactionRates());

		writeMapAsColumnFile(prefix + "heat.dat", gd.otherHeating());
		writeMapAsColumnFile(prefix + "cool.dat", gd.cooling());

		ColumnFile radfield(prefix + "nu_jnu.dat", {"frequency", "nu Jnu"});
		for (size_t i = 0; i < frequencyv.size(); i++)
		{
			double wav = Constant::LIGHT / frequencyv[i];
			radfield.writeLine<Array>({wav * Constant::CM_UM,
			                           Constant::FPI * frequencyv[i] *
			                                           specificIntensityv[i]});
		}
		GasModule::GasState gs = s.makeGasState(gasInterface.oFrequencyv(),
		                                        gasInterface.eFrequencyv());
		writeGasState(prefix, gasInterface, gs);
		writeGrains(prefix, gri);
		// plotHeatingCurve(*gasInterface.pimpl(), "MRNDust/", nHtotal, I_nu, gri);

		vector<string> populationColumns = {"label", "energy", "density"};
		ColumnFile hpop(prefix + "hpopulations.dat", populationColumns);
		for (double n : gd.hPopulations())
			hpop.writeLine<Array>({0., 0., n});
		ColumnFile h2pop(prefix + "h2populations.dat", populationColumns);
		for (double n : gd.h2Populations())
			h2pop.writeLine<Array>({0., 0., n});
	}
}
