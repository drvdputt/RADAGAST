#include "H2FromFiles.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"

#include <regex>

using namespace std;

namespace
{
constexpr bool BLEVELS{true};
constexpr bool CLEVELS{true};
} // namespace

H2FromFiles::H2FromFiles(int maxJ, int maxV)
                : _maxJ{maxJ}, _maxV{maxV}, _inH{SpeciesIndex::inH()},
                  _inH2{SpeciesIndex::inH2()}
{
	readLevels();
	readTransProbs();
	readCollisions();
	readDissProbs();
	readDirectDissociation();
// #define PLOT_LEVEL_MATRICES
#ifdef PLOT_LEVEL_MATRICES
	ofstream avvOut = IOTools::ofstreamFile("h2/einsteinA.dat");
	avvOut << _avv << endl;
	avvOut.close();

	ofstream lvlOut = IOTools::ofstreamFile("h2/levels.dat");
	lvlOut << "# eState\tv\tJ\tE\n";
	for (const auto& l : _levelv)
		lvlOut << static_cast<int>(l.eState()) << "\t" << l.v() << "\t" << l.j() << "\t"
		       << l.e() << endl;
	lvlOut.close();
#endif
}

size_t H2FromFiles::numLv() const { return _numL; }

size_t H2FromFiles::indexOutput(ElectronicState eState, int j, int v) const
{
	// This cast is safe, as the default underlying type of an enum is int.
	return _ejvToIndexm.at({static_cast<int>(eState), j, v});
}

void H2FromFiles::readLevels()
{
	// Expand levelv with the levels listed in these files
	readLevelFile("dat/h2/energy_X.dat", ElectronicState::X);
	if (BLEVELS)
		readLevelFile("dat/h2/energy_B.dat", ElectronicState::B);
	if (CLEVELS)
	{
		readLevelFile("dat/h2/energy_C_plus.dat", ElectronicState::Cplus);
		readLevelFile("dat/h2/energy_C_minus.dat", ElectronicState::Cminus);
	}
	_numL = _levelv.size();
}

void H2FromFiles::readTransProbs()
{
	_avv = EMatrix::Zero(_numL, _numL);

	// Radiative transitions within ground stated. Fills _avv with data from Wolniewicz
	// (1998)
	readTransProbFile("dat/h2/transprob_X.dat", ElectronicState::X, ElectronicState::X);
	if (BLEVELS)
		readTransProbFile("dat/h2/transprob_B.dat", ElectronicState::B,
		                  ElectronicState::X);
	if (CLEVELS)
	{
		readTransProbFile("dat/h2/transprob_C_plus.dat", ElectronicState::Cplus,
		                  ElectronicState::X);
		readTransProbFile("dat/h2/transprob_C_minus.dat", ElectronicState::Cminus,
		                  ElectronicState::X);
	}
}

void H2FromFiles::readDissProbs()
{
	_dissProbv = EVector::Zero(_numL);
	_dissKinEv = EVector::Zero(_numL);

	if (BLEVELS)
		readDissProbFile("dat/h2/dissprob_B.dat", ElectronicState::B);
	if (CLEVELS)
	{
		readDissProbFile("dat/h2/dissprob_C_plus.dat", ElectronicState::Cplus);
		readDissProbFile("dat/h2/dissprob_C_minus.dat", ElectronicState::Cminus);
	}
}

void H2FromFiles::readCollisions()
{
	// Collisions with H, from Lique 2015
	_qH = readCollisionFile("dat/h2/coll_rates_H_15.dat");

	// Collision rates with H2, from Lee 2008
	_qH2ortho = readCollisionFile("dat/h2/coll_rates_H2ortho_ORNL.dat");
	_qH2para = readCollisionFile("dat/h2/coll_rates_H2para_ORNL.dat");
}

void H2FromFiles::readDirectDissociation()
{
	_dissociationCrossSectionv.resize(_levelv.size());
	ifstream cont_diss = IOTools::ifstreamRepoFile("dat/h2/cont_diss.dat");
	string line;

	// Main file traversal loop
	while (getline(cont_diss, line))
	{
		// Find the start of a cross section data block
		if (line.at(0) == '#' && line.at(1) == '!')
		{
			// Get the three lines starting with "#!". I don't think I need the
			// 2nd and 3rd lines though.
			string line1, line2;
			getline(cont_diss, line1);
			getline(cont_diss, line2);

			// Parse the first line to get the quantum numbers
			int nei, nef, vi, ji;
			string cleanline = regex_replace(line, regex("[^0-9.]+"), " ");
			istringstream(cleanline) >> nei >> nef >> vi >> ji;

			// Skip to the next #! block if we are not treating this level (There
			// should be no data for states other than X, but let's ignore nei
			// anyway. .
			if (nei > 0 || !validJV(ji, vi))
				continue; // This will never be used;

			// Read the data
			string dataline;
			vector<double> frequencyv, crossSectionv;

			while (getline(cont_diss, dataline))
			{
				// Parse the data line:	Split by comma
				istringstream iss(dataline);
				string s0, s1;
				getline(iss, s0, ',');
				getline(iss, s1);

				// String to double
				double energy_invcm = stod(s0);
				double crossSection_ang2 = stod(s1);

				// Convert to the correct units and add
				frequencyv.emplace_back(Constant::LIGHT * energy_invcm);
				crossSectionv.emplace_back(crossSection_ang2 *
				                           Constant::ANG_CM * Constant::ANG_CM);

				if (cont_diss.peek() == '#')
					break;
			}
			// When this loop exits, we should be at the start of the next comment
			// block, or at the end of the file.

			// Add the new cross section at the correct level index.
			size_t levelIndex = indexOutput(ElectronicState::X, ji, vi);
			_dissociationCrossSectionv[levelIndex].emplace_back(
			                Array(frequencyv.data(), frequencyv.size()),
			                Array(crossSectionv.data(), crossSectionv.size()));
		}
	}
}

void H2FromFiles::readLevelFile(const string& repoFile, ElectronicState eState)
{
	ifstream energy = IOTools::ifstreamRepoFile(repoFile);

	string line;
	getline(energy, line);
	int y, m, d;
	istringstream(line) >> y >> m >> d;
	Error::equalCheck("magic number y", y, 2);
	Error::equalCheck("magic number m", m, 4);
	Error::equalCheck("magic number d", d, 29);

	// Start reading in the V, J and E (cm-1) of the levels
	vector<int> vv, jv;
	vector<double> ev;
	int counter = 0;
	while (getline(energy, line))
	{
		auto iss = istringstream(line);
		if (line.front() == '#')
		{
			string word1, word2;
			iss >> word1 >> word2;
			// If the comment is anything else than "#extra data", we will skip this
			// line
			if (word1 != "#extra" || word2 != "data")
				continue;
		}
		int v, j;
		iss >> v >> j;
		double k; // cm-1
		iss >> k;

		// Add the level if within the requested J,v limits. Convert 1/wavelength [cm-1]
		// to ergs by using E = hc / lambda.
		if (validJV(j, v))
		{
			counter++;
			addLevel(eState, j, v, Constant::PLANCKLIGHT * k);
		}
	}
	DEBUG("Read in " << counter << " levels from " << repoFile << endl);
}

void H2FromFiles::readTransProbFile(const string& repoFile, ElectronicState upperE,
                                    ElectronicState lowerE)
{
	ifstream transprob = IOTools::ifstreamRepoFile(repoFile);

	string line;
	getline(transprob, line);
	int y, m, d;
	istringstream(line) >> y >> m >> d;
	Error::equalCheck("magic number y", y, 2);
	Error::equalCheck("magic number m", m, 4);
	Error::equalCheck("magic number d", d, 29);

	size_t counter = 0;
	while (getline(transprob, line))
	{
		if (line.empty() || line.at(0) == '#')
			continue;

		int EU, VU, JU, EL, VL, JL;
		double A; // Unit s-1
		istringstream(line) >> EU >> VU >> JU >> EL >> VL >> JL >> A;

		// Check if the entry in the file and the expected Estate-changing transition
		// correspond (just a double check on the enum scheme)
		Error::equalCheck("EU and cast of upperE", EU, static_cast<int>(upperE));
		Error::equalCheck("EL and cast of lowerE", EL, static_cast<int>(lowerE));

		// If within the J,v limits, fill in the coefficient.
		if (validJV(JU, VU) && validJV(JL, VL))
		{
			// An error will be thrown if the level is still not in the list.
			size_t upperIndex = indexOutput(upperE, JU, VU);
			size_t lowerIndex = indexOutput(lowerE, JL, VL);
			_avv(upperIndex, lowerIndex) = A;
			counter++;
		}
	}
	DEBUG("Read in " << counter << " Einstein A coefficients from " << repoFile << endl);
}

void H2FromFiles::readDissProbFile(const string& repoFile, ElectronicState eState)
{
	ifstream dissprobs = IOTools::ifstreamRepoFile(repoFile);
	int y, m, d;
	IOTools::istringstreamNextLine(dissprobs) >> y >> m >> d;
	Error::equalCheck("magic y", y, 3);
	Error::equalCheck("magic m", m, 2);
	Error::equalCheck("magic d", d, 11);

	string line;
	int counter = 0;
	while (getline(dissprobs, line))
	{
		if (line.empty() || line.at(0) == '#')
			continue;

		int v, j;
		// diss prob in s-1, kin energy in eV
		double diss, kin;
		istringstream(line) >> v >> j >> diss >> kin;

		if (!validJV(j, v))
			continue;

		int index = indexOutput(eState, j, v);
		_dissProbv[index] = diss;
		_dissKinEv[index] = kin / Constant::ERG_EV;
		counter++;
	}
	DEBUG("Read in " << counter << " dissociation rates from " << repoFile << endl);
}

CollisionData H2FromFiles::readCollisionFile(const string& repoFile) const
{
	ifstream coll_rates = IOTools::ifstreamRepoFile(repoFile);

	int m;
	IOTools::istringstreamNextLine(coll_rates) >> m;
	Error::equalCheck(repoFile + "magic number", m, 110416);

	// Read until the temperature line
	string line;
	while (getline(coll_rates, line))
		if (line.at(0) != '#')
			break;

	// The number of dots gives us the number of columns
	size_t numTemperatures = count(begin(line), end(line), '.');

	// Read in the temperatures
	Array temperaturev(numTemperatures);
	istringstream issTemperatures(line);
	for (size_t i = 0; i < numTemperatures; i++)
		issTemperatures >> temperaturev[i];

	// For each line, get the initial level index, final level index, and the data
	vector<int> iv, fv;
	vector<Array> qvv;
	while (getline(coll_rates, line))
	{
		if (line.empty() || line.at(0) == '#')
			continue;

		istringstream issCollRates(line);
		int VU, JU, VL, JL;
		issCollRates >> VU >> JU >> VL >> JL;

		// If we are not treating this level (the J or the v is out of range), skip.
		if (!validJV(JU, VU) || !validJV(JL, VL))
			continue;

		// Else, add the level indices and data.
		iv.emplace_back(indexOutput(ElectronicState::X, JU, VU));
		fv.emplace_back(indexOutput(ElectronicState::X, JL, VL));
		Array qForEachTv(numTemperatures);
		for (size_t i = 0; i < numTemperatures; i++)
			issCollRates >> qForEachTv[i];
		qvv.emplace_back(qForEachTv);
	}

	size_t numTransitions = qvv.size();
	CollisionData qData;
	qData.prepare(temperaturev, numTransitions);
	for (size_t t = 0; t < numTransitions; t++)
		qData.insertDataForTransition(qvv[t], iv[t], fv[t]);
	qData.check();
	DEBUG("Read in " << qData.transitionv().size() << " collision coefficients from "
	                 << repoFile << endl);
	return qData;
}

EVector H2FromFiles::ev() const
{
	EVector the_ev(_numL);
	for (size_t i = 0; i < _numL; i++)
		the_ev(i) = _levelv[i].e();
	return the_ev;
}

EVector H2FromFiles::gv() const
{
	EVector the_gv(_numL);
	for (size_t i = 0; i < _numL; i++)
		the_gv(i) = _levelv[i].g();
	return the_gv;
}

EMatrix H2FromFiles::avv() const { return _avv; }

EMatrix H2FromFiles::extraAvv() const { return EMatrix::Zero(_numL, _numL); }

EMatrix H2FromFiles::cvv(double T, const EVector& speciesNv) const
{
	EMatrix the_cvv{EMatrix::Zero(_numL, _numL)};

	// H-H2 collisions
	double nH = speciesNv(_inH);
	addToCvv(the_cvv, _qH, T, nH);

	// TODO: (optional?) g-bar approximation for missing coefficients

	// H2-H2 collisions. TODO: I actually need the ortho-para ratio here, which
	// unfortunately depends on the levels themselves... Maybe I need a wrapper around
	// speciesNv, which also contains a bunch of other parameters, such as the ortho and
	// para contributions. I just pick a ratio of 3 here, hardcoded.
	double nH2 = speciesNv(_inH2);
	addToCvv(the_cvv, _qH2ortho, T, 0.75 * nH2);
	addToCvv(the_cvv, _qH2para, T, 0.25 * nH2);

	return the_cvv;
}

double H2FromFiles::directDissociationCrossSection(double nu, int j, int v) const
{
	return directDissociationCrossSection(nu, indexOutput(ElectronicState::X, j, v));
}

double H2FromFiles::directDissociationCrossSection(double nu, size_t index) const
{
	double sigma{0.};
	// Evaluate all the cross sections for this level at this frequency
	for (const Spectrum& cs : _dissociationCrossSectionv[index])
		sigma += cs.evaluate(nu);
	return sigma;
}

vector<Spectrum> H2FromFiles::directDissociationCrossSections(size_t index) const
{
	return _dissociationCrossSectionv[index];
}

bool H2FromFiles::validJV(int J, int v) const { return J <= _maxJ && v <= _maxV; }

void H2FromFiles::addToCvv(EMatrix& the_cvv, const CollisionData& qdata, double T,
                           double nPartner) const
{
	if (nPartner <= 0)
		return;

	/* Find the grid point to the right of (>=) the requested log-temperature. (Returns last
	   point if T > Tmax). We will naively extrapolate for points outside the range, and cut
	   off the result at 0 if it becomes negative. */
	const Array& temperaturev = qdata.temperaturev();
	int iRight = TemplatedUtils::index(T, temperaturev);
	iRight = max(iRight, 1);
	int iLeft = iRight - 1;
	double tRight = temperaturev[iRight];
	double tLeft = temperaturev[iLeft];

	double kT = T * Constant::BOLTZMAN;
	vector<array<int, 2>> transitionv = qdata.transitionv();
	for (int transitionIndex = 0; transitionIndex < transitionv.size(); transitionIndex++)
	{
		// Interpolate data points naively for temperatures left and right of T
		double qLeft = qdata.q(iLeft, transitionIndex);
		double qRight = qdata.q(iRight, transitionIndex);
		double q = TemplatedUtils::interpolateLinear(T, tLeft, tRight, qLeft, qRight);
		// Make zero if interpolation makes this negative
		q = max(0., q);

		// The levels involved in this transition
		int i = transitionv[transitionIndex][0];
		int f = transitionv[transitionIndex][1];
		const H2Level& ini = _levelv[i];
		const H2Level& fin = _levelv[f];
		double Cif = q * nPartner;
		double Cfi = Cif * ini.g() / fin.g() * exp((fin.e() - ini.e()) / kT);
		the_cvv(i, f) += Cif;
		the_cvv(f, i) += Cfi;
	}
}
