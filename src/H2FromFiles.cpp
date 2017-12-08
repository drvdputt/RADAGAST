#include "H2FromFiles.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"

// TODO: make model

using namespace std;

H2FromFiles::H2FromFiles(int maxJ, int maxV)
                : _maxJ{maxJ}, _maxV{maxV}, _inH{SpeciesIndex::inH()}
{
	readLevels();
	readTransProb();
	readCollisions();
#ifdef PRINT_LEVEL_MATRICES
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
	ifstream energyX = IOTools::ifstreamRepoFile("dat/h2/energy_X.dat");
	string line;

	// Get the date / magic number
	getline(energyX, line);
	int y, m, d;
	istringstream(line) >> y >> m >> d;
	Error::equalCheck("y", y, 2);
	Error::equalCheck("m", m, 4);
	Error::equalCheck("d", d, 29);

	// Skip the next two lines;
	getline(energyX, line);
	getline(energyX, line);

	// The eState will be fixed to X here
	auto eState{ElectronicState::X};
	// Start reading in the V, J and E (cm-1) of the levels
	_levelv.reserve(301);
	while (getline(energyX, line))
	{
		auto iss = istringstream(line);

		// If there is a comment in front of the entry, extract it as a string.
		if (line.front() == '#')
		{
			string word1, word2;
			iss >> word1 >> word2;
			// If the comment is anything else than "#extra data", we will skip this
			// line
			if (word1 != "#extra" || word2 != "data")
				continue;
		}

		// Get the numbers
		int v, j;
		iss >> v >> j;
		// Wave number (in cm-1)
		double k;
		iss >> k;

		// Add the level if within the given limits.
		if (validJV(j, v))
			addLevel(eState, j, v, Constant::PLANCKLIGHT * k);
	}
	_numL = _levelv.size();
	DEBUG("Read in " << _numL << " H2 levels" << endl);
}

void H2FromFiles::readTransProb()
{
	_avv = EMatrix::Zero(_numL, _numL);

	ifstream transprobX = IOTools::ifstreamRepoFile("dat/h2/transprob_X.dat");
	string line;
	getline(transprobX, line);
	int y, m, d;
	istringstream(line) >> y >> m >> d;

	// Transition for eState X (ground state) (Data from Wolniewicz (1998))
	auto eState{ElectronicState::X};
	size_t counter{0};
	while (getline(transprobX, line))
	{
		if (line.empty() || line.at(0) == '#')
			continue;

		int EU, VU, JU, EL, VL, JL;
		double A; // Unit s-1
		istringstream(line) >> EU >> VU >> JU >> EL >> VL >> JL >> A;
		/* TODO: check if the energy of level upperIndex is actually higher than energy of
		   level lowerIndex */
		// Dont use EU and EL here, as we already know that it's all for X

		// If within the limits, fill in the coefficient.
		if (validJV(JU, VU) && validJV(JL, VL))
		{
			// An error will be thrown if the level is still not in the list.
			size_t upperIndex = indexOutput(eState, JU, VU);
			size_t lowerIndex = indexOutput(eState, JL, VL);
			_avv(upperIndex, lowerIndex) = A;
			counter++;
		}
	}
	DEBUG("Read in " << counter << " Einstein A coefficients." << endl);
}

void H2FromFiles::readCollisions()
{
	ifstream coll_rates_H_15 = IOTools::ifstreamRepoFile("dat/h2/coll_rates_H_15.dat");

	int m;
	IOTools::istringstreamNextLine(coll_rates_H_15) >> m;
	Error::equalCheck("coll_rates_H_15 magic numbers", m, 110416);

	// Read until the temperature line
	string line;
	while (getline(coll_rates_H_15, line))
		if (line.at(0) != '#')
			break;

	// Hardcode this
	const size_t numTemperatures{50};
	Array temperaturev(numTemperatures);
	istringstream tempLine(line);
	for (size_t i = 0; i < numTemperatures; i++)
		tempLine >> temperaturev[i];

	// Do not care about efficiency for now.
	// For each line, get the initial level index, final level index, and the data.
	vector<int> iv, fv;
	vector<Array> qvv; // indexed on [transition][temperature]
	while (getline(coll_rates_H_15, line))
	{
		istringstream transitionLine(line);
		int VU, JU, VL, JL;
		transitionLine >> VU >> JU >> VL >> JL;

		// If we are not treating one of these (J,v)'s, skip to the next line
		if (!validJV(JU, VU) || !validJV(JL, VL))
			continue;

		// Else, add the level indices and data.
		iv.emplace_back(indexOutput(ElectronicState::X, JU, VU));
		fv.emplace_back(indexOutput(ElectronicState::X, JL, VL));
		Array qForEachTv(numTemperatures);
		for (size_t i = 0; i < numTemperatures; i++)
			transitionLine >> qForEachTv[i];
		qvv.emplace_back(qForEachTv);
	}

	// Now put this information into the CollisionData object
	size_t numTransitions = qvv.size();
	_qH.prepare(temperaturev, numTransitions);
	for (size_t t = 0; t < numTransitions; t++)
		_qH.insertDataForTransition(qvv[t], iv[t], fv[t]);
	_qH.check();
	DEBUG("Read in " << _qH.transitionv().size() << " collision coefficients." << endl);
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

	return the_cvv;
}

EVector H2FromFiles::sourcev(double T, const EVector& speciesNv) const
{
	return EVector::Zero(_numL);
}

EVector H2FromFiles::sinkv(double T, double n, const EVector& speciesNv) const
{
	return EVector::Zero(_numL);
}

bool H2FromFiles::validJV(int J, int v) const { return J <= _maxJ && v <= _maxV; }

void H2FromFiles::addToCvv(EMatrix& the_cvv, const CollisionData& qdata, double T,
                           double nPartner) const
{
	if (nPartner <= 0)
		return;

	/* Find the grid point to the right of (>=) the requested log-temperature. (Returns last
	   point if T > Tmax). We will naively extrapolate for points outside the range, and cut off
	   the result at 0 if it becomes negative. */
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
