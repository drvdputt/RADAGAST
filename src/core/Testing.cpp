#include "Testing.h"
#include "ChemicalNetwork.h"
#include "ChemistrySolver.h"
#include "FreeBound.h"
#include "GasInterface.h"
#include "GasInterfaceImpl.h"
#include "GasStruct.h"
#include "GrainType.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenHardcoded.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SpecialFunctions.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"

#include <iomanip>
#include <iterator>
#include <sstream>

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
	size_t aIndex = TemplatedUtils::index(a, qds.FILEAV);

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

#define PLOT_QABS
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

Array Testing::improveFrequencyGrid(const NLevel& boundBound, const Array& oldPoints)
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

Array Testing::defaultCoarseFrequencyv()
{
	Array coarsev = generateGeometricGridv(100, defaultMinFreq, defaultMaxFreq);
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
                                  const std::string& outputPath, double Tc, double G0, double n,
                                  double expectedTemperature)
{
	Array specificIntensityv = generateSpecificIntensityv(gi.iFrequencyv(), Tc, G0);

	GasModule::GasState gs;
	GasModule::GrainInterface gri{};
	gi.updateGasState(gs, n, expectedTemperature, specificIntensityv, gri);
	writeGasState(outputPath, gi, gs);
}

void Testing::writeGasState(const string& outputPath, const GasModule::GasInterface& gi,
                            const GasModule::GasState& gs)
{
	cout << "Equilibrium temperature: " << gs.temperature() << endl;
	cout << "Ionized fraction: " << gs.ionizedFraction() << endl;

	const Array& eFrequencyv = gi.eFrequencyv();
	const Array& emv = gs._emissivityv;
	Spectrum opv(gi.oFrequencyv(), gs._opacityv);

	cout << "Integrated emissivity " << TemplatedUtils::integrate<double>(eFrequencyv, emv)
	     << endl;

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
		out << scientific << freq << tab << wav << tab << emv[iFreq] << tab
		    << opv.evaluate(freq) << endl;
		wavfile.precision(9);
		wavfile << wav << tab << freq << endl;
	}
	out.close();
	wavfile.close();

	out = IOTools::ofstreamFile(outputPath + "raw_opacity.dat");
	for (size_t i = 0; i < gs._opacityv.size(); i++)
	{
		out << gi.oFrequencyv()[i] << tab << gs._opacityv[i] << endl;
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
                               const GasModule::GrainInterface& gri)
{
	const string tab = "\t";

	ofstream output = IOTools::ofstreamFile(outputPath + "heatingcurve.dat");
	output << "# 0temperature 1net 2heat 3cool"
	       << " 4lineNet 5lineHeat 6lineCool"
	       << " 7continuumNet 8continuumHeat 9continuumCool"
	       << "10e- 11H+ 12H 13H2 \n";

	double T = 10;
	Array Tv = Testing::generateGeometricGridv(100, 10, 50000);
	GasInterfaceImpl::Solution s =
	                gi.calculateDensities(n, T, specificIntensity, gri, nullptr);

	auto outputCooling = [&](double t) {
		s = gi.calculateDensities(n, t, specificIntensity, gri, &s);
		double heat = gi.heating(s);
		double cool = gi.cooling(s);
		double lHeat = gi.lineHeating(s);
		double lCool = gi.lineCooling(s);
		double cHeat = gi.continuumHeating(s);
		double cCool = gi.continuumCooling(s);

		double netHeating = heat - cool;
		double netLine = lHeat - lCool;
		double netCont = cHeat - cCool;

		output << t << tab << netHeating << tab << heat << tab << cool << tab << netLine
		       << tab << lHeat << tab << lCool << tab << netCont << tab << cHeat << tab
		       << cCool;
		for (int i = 0; i < s.speciesNv.size(); i++)
			output << tab << s.speciesNv(i);
		output << '\n';
	};

	// forward sweep
	for (double& T : Tv)
		outputCooling(T);
	// reverse sweep
	size_t i = Tv.size();
	while (i--)
		outputCooling(Tv[i]);

	output.close();

	double isrf = TemplatedUtils::integrate<double>(specificIntensity.frequencyv(),
	                                                specificIntensity.valuev());

	cout << "Calculated heating curve under isrf of " << isrf << " erg / s / cm2 / sr = "
	     << isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl;
}

void Testing::plotPhotoelectricHeating()
{
	bool carbon{true};
	readQabs(carbon);

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

	HydrogenFromFiles hff(5);
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
			size_t nliIndex = hff.indexOutput(n, li);
			// Decay rate to all other levels
			double anl = anlv(nliIndex);
			double qnl = 0;

			// Get the total decay rate due to l-changing collisions
			for (int lf = 0; lf < n; lf++)
			{
				size_t nlfIndex = hff.indexOutput(n, lf);
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
	H2Levels h2l(make_shared<H2FromFiles>(maxJ, maxV));
	Array frequencyv = improveFrequencyGrid(h2l, unrefinedv);

	auto h2Data{make_shared<H2FromFiles>(maxJ, maxV)};
	H2Levels h2Levels{h2Data};

	// Set the densities
	EVector speciesNv{EVector::Zero(SpeciesIndex::size())};
	speciesNv(SpeciesIndex::inH2()) = nH2;
	speciesNv(SpeciesIndex::ine()) = ne;
	speciesNv(SpeciesIndex::inp()) = np;
	speciesNv(SpeciesIndex::inH()) = nH;

	GasStruct gas(T, speciesNv);
	NLevel::Solution s = h2l.customSolution(nH2, gas, specificIntensity);
	if (write)
	{
		Array emissivityv = h2Levels.emissivityv(s, frequencyv);
		Array opacityv = h2Levels.opacityv(s, unrefinedv);
		Array lineOp = h2Levels.lineOpacityv(s, unrefinedv);

		ofstream h2optical = IOTools::ofstreamFile("h2/opticalProperties.dat");
		h2optical << "# nu (hz)\tlambda (micron)\temissivity\topacity\tlineOpacity\n";
		const string tab{"\t"};
		for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		{
			h2optical << frequencyv[iFreq] << tab
			          << Constant::LIGHT / frequencyv[iFreq] * Constant::CM_UM
			          << tab << emissivityv[iFreq] << tab << opacityv[iFreq] << tab
			          << lineOp[iFreq] << endl;
		}
	}
}

void Testing::runFromFilesvsHardCoded()
{
	Array unrefinedv =
	                generateGeometricGridv(1000, Constant::LIGHT / (1e10 * Constant::UM_CM),
	                                       Constant::LIGHT / (0.00001 * Constant::UM_CM));

	// Hey, at least we'll get a decent frequency grid out of this hack
	HydrogenLevels hl(make_shared<HydrogenFromFiles>(5));
	FreeBound fb;
	Array frequencyv = improveFrequencyGrid(hl, unrefinedv);
	frequencyv = improveFrequencyGrid(fb, frequencyv);

	GasModule::GasInterface gihhc(unrefinedv, unrefinedv, frequencyv, "hhc", "none");
	runGasInterfaceImpl(gihhc, "hardcoded/");

	GasModule::GasInterface gihff(unrefinedv, unrefinedv, frequencyv, "hff2", "none");
	runGasInterfaceImpl(gihff, "fromfiles/");
}

GasModule::GasInterface Testing::genFullModel()
{
	Array coarsev = defaultCoarseFrequencyv();

	cout << "Construction model to help with refining frequency grid" << endl;
	HydrogenLevels hl(make_shared<HydrogenFromFiles>());
	FreeBound fb;
	H2Levels h2l(make_shared<H2FromFiles>(8, 5));

	Array frequencyv = improveFrequencyGrid(hl, coarsev);
	frequencyv = improveFrequencyGrid(fb, frequencyv);
	frequencyv = improveFrequencyGrid(h2l, frequencyv);

	cout << "Constructing new model using the improved frequency grid" << endl;
	return {coarsev, coarsev, frequencyv, "", "8 5"};
}

GasModule::GasInterface Testing::genHonlyModel()
{
	Array coarsev = defaultCoarseFrequencyv();

	cout << "Construction model to help with refining frequency grid" << endl;
	HydrogenLevels hl(make_shared<HydrogenFromFiles>());
	FreeBound fb;

	Array frequencyv = improveFrequencyGrid(hl, coarsev);
	frequencyv = improveFrequencyGrid(fb, frequencyv);

	cout << "Constructing new model using the improved frequency grid" << endl;
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

void Testing::runWithDust(bool write)
{
	cout << "RUN_WITH_DUST\n";

	// Gas model
	GasModule::GasInterface gasInterface = genFullModel();

	// Radiation field
	double Tc{4e3};
	double G0{1e2};
	Array specificIntensityv =
	                generateSpecificIntensityv(gasInterface.iFrequencyv(), Tc, G0);

	// Gas density
	double nHtotal{1000};
	double Tinit{8000};

	// TODO: need a reasonable grain size distribution and DGR for testing here
	Array sizev, densityv, temperaturev;
	// Provide sizes in cm
	// 10 A, 100A ,1000 A
	sizev = {1e-7, 1e-6, 1e-5};
	densityv = {1e-5, 1e-8, 1e-11};
	// number of grains per H atom (this should be very small, as grains are much heavier
	// than a hydrogen atom AND their total mass is only about 1% of the hydrogen mass

	densityv *= nHtotal; // density in cm-3
	// densityv *= 0;

	// Actually, a temperature distribution per grain size is the most detailed form we'll
	// be using. Maybe we need multiple versions of the grain interface with regards to
	// storing the grain temperatures.
	temperaturev = {15, 10, 5};
	// And now I need some absorption efficiencies for every wavelength. Let's try to use
	// the old photoelectric heating test code.
	bool car = true;
	vector<Array> qAbsvv{qAbsvvForTesting(car, sizev, gasInterface.iFrequencyv())};

	// TODO: check if the qabsvv has loaded correctly

	// Construct grain info using list of population objects
	auto grainPopv{make_unique<vector<GasModule::GrainInterface::Population>>()};
	grainPopv->emplace_back(GasModule::GrainTypeLabel::CAR, sizev, densityv, temperaturev,
	                        qAbsvv);
	GasModule::GrainInterface grainInterface(move(grainPopv));

	// Run
	GasModule::GasState gs;
	gasInterface.updateGasState(gs, nHtotal, Tinit, specificIntensityv, grainInterface);
	if (write)
		writeGasState("withDust/", gasInterface, gs);
}

// void Testing::plotHlines()
// {
//	HydrogenFromFiles hdata{};
//	HydrogenLevels hlevels{make_shared<HydrogenDataProvider>(hdata)};
//	double Tc = 30000;
//	double g0 = 1e0;
//	double n = 1000;
//	Array frequencyv = defaultCoarseFrequencyv();
//	Spectrum specificIntensity{frequencyv, generateSpecificIntensityv(frequencyv, Tc, g0)};
//	NLevel::Solution s;
//	s.n = n;
//	s.nv = EVector::Zero(SpeciesIndex::size());
//	s.nv(SpeciesIndex::inH()) = n;
//	EMatrix Tvv = hlevels.totalTransitionRatesvv(specificIntensity,const GasStruct& gas,EMatrix* cvv_p[=nullptr])
//	Array lineSpectrum = hlevels.lineEmissivityv(s, frequencyv);
// }

GasModule::GrainInterface Testing::genMRNDust(double nHtotal, const Array& frequencyv)
{
	auto grainPopv{make_unique<vector<GasModule::GrainInterface::Population>>()};

	double amin = 50 * Constant::ANG_CM;
	double amax = 0.24 * Constant::UM_CM;

	// Shared properties
	size_t numSizes = 100;
	Array bin_edges = generateGeometricGridv(numSizes + 1, amin, amax);
	Array sizev(numSizes);
	Array temperaturev(numSizes);
	for (size_t i = 0; i < numSizes; i++)
	{
		sizev[i] = (bin_edges[i + 1] + bin_edges[i]) / 2;
		temperaturev[i] = 15;
	}

	// Properties that differ between car and sil
	auto addMRNCarOrSil = [&](bool car) {
		double log10norm = car ? -25.13 : -25.11;
		double C = pow(10., log10norm);
		Array densityv(numSizes);
		for (size_t i = 0; i < numSizes; i++)
		{
			double bin_width = bin_edges[i + 1] - bin_edges[i];
			densityv[i] = C * nHtotal * pow(sizev[i], -3.5) * bin_width;
		}
		auto label = car ? GasModule::GrainTypeLabel::CAR
		                 : GasModule::GrainTypeLabel::SIL;
		const auto& qabsvv = qAbsvvForTesting(car, sizev, frequencyv);
		grainPopv->emplace_back(label, sizev, densityv, temperaturev, qabsvv);
	};

	// Carbonaceous population:
	addMRNCarOrSil(true);

	// Silicate population:
	addMRNCarOrSil(false);

	GasModule::GrainInterface grainInterface(move(grainPopv));
	return grainInterface;
}

void Testing::runMRNDust(bool write)
{
	cout << "RUN_MRN_DUST\n";

	// Gas model
	GasModule::GasInterface gasInterface = genFullModel();
	double nHtotal{1000};
	double Tinit{8000};

	// Dust model
	auto gri = genMRNDust(nHtotal, gasInterface.iFrequencyv());

	// Radiation field
	double Tc{4e3};
	double G0{1e2};
	Array specificIntensityv =
	                generateSpecificIntensityv(gasInterface.iFrequencyv(), Tc, G0);

	GasModule::GasState gs;
	gasInterface.updateGasState(gs, nHtotal, Tinit, specificIntensityv, gri);

	if (write)
	{
		writeGasState("MRNDust/", gasInterface, gs);
		Spectrum I_nu = Spectrum(gasInterface.iFrequencyv(), specificIntensityv);
		plotHeatingCurve(*gasInterface.pimpl(), "MRNDust/", nHtotal, I_nu, gri);
	}
}
