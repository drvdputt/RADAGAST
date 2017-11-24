#include "H2FromFiles.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"

// TODO: make model

using namespace std;

const string dataLocation{REPOROOT "/dat/h2"};

H2FromFiles::H2FromFiles() { readData(); }

size_t H2FromFiles::numLv() const { return _numL; }

size_t H2FromFiles::indexOutput(ElectronicState eState, int j, int v) const
{
	// This cast is safe, as the default underlying type of an enum is int.
	return _evjToIndex.at({static_cast<int>(eState), j, v});
}

void H2FromFiles::readData()
{
	//-----------------//
	// READ LEVEL DATA //
	//-----------------//
	ifstream energyX = IOTools::ifstreamFile(dataLocation + "/energy_X.dat");
	string line;

	// Get the date / magic number
	getline(energyX, line);
	int y, m, d;
	istringstream(line) >> y >> m >> d;
	Error::equalCheck("y", y, 2);
	Error::equalCheck("m", m, 4);
	Error::equalCheck("d", d, 29);

	// Skip the next two lines
	IOTools::skipLines(energyX, 2);

	// The eState will be fixed to X here
	ElectronicState eState{ElectronicState::X};
	// Start reading in the V, J and E (cm-1) of the levels
	_levelv.reserve(301);
	while (getline(energyX, line))
	{
		auto iss = istringstream(line);

		// If there is a comment in front of the entry, extract it as a string
		if (line.front() == '#')
		{
			string word1, word2;
			iss >> word1 >> word2;
			DEBUG(word1 + word2 << endl);
		}

		// Get the numbers
		int v, j;
		iss >> v >> j;
		// Wave number
		double k;
		iss >> k;

		addLevel(eState, j, v, Constant::PLANCKLIGHT * k);
		DEBUG(_levelv.size() << " " << j << " " << v << " "
		                     << _levelv.back().e() * Constant::ERG_EV << " eV" << endl);
	}
	_numL = _levelv.size();
	DEBUG("Read in " << _numL << " H2 levels" << endl);
	energyX.close();

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
	_avv = EMatrix::Zero(_numL, _numL);

	ifstream transprobX = IOTools::ifstreamFile(dataLocation + "/transprob_X.dat");
	getline(transprobX, line);
	istringstream(line) >> y >> m >> d;

	// Transition for eState X (ground state) (Data from Wolniewicz (1998))
	eState = ElectronicState::X;
	while (getline(transprobX, line))
	{
		// Skip comments or empty lines
		if (line.empty() || line.at(0) == '#')
			continue;

		int EU, VU, JU, EL, VL, JL;
		double A; // Unit s-1
		istringstream(line) >> EU >> VU >> JU >> EL >> VL >> JL >> A;
		/* TODO: check if the energy of level upperIndex is actually higher than energy of
		   level lowerIndex */
		// Dont use EU and EL here, as we already know that it's all for X
		size_t upperIndex = indexOutput(eState, JU, VU);
		size_t lowerIndex = indexOutput(eState, JL, VL);
		_avv(upperIndex, lowerIndex) = A;
	}
	transprobX.close();
#ifdef PRINT_LEVEL_MATRICES
	ofstream avvOut = IOTools::ofstreamFile("h2/einsteinA.dat");
	avvOut << _avv << endl;
	avvOut.close();
#endif
#define PLOT_ENERGIES
#ifdef PLOT_ENERGIES
	ofstream lvlOut = IOTools::ofstreamFile("h2/levels.dat");
	lvlOut << "# eState\tv\tJ\tE\n";
	for (const auto& l : _levelv)
		lvlOut << static_cast<int>(l.eState()) << "\t" << l.v() << "\t" << l.j() << "\t"
		       << l.e() << endl;
	lvlOut.close();
#endif
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
	return EMatrix::Zero(_numL, _numL);
}

EVector H2FromFiles::sourcev(double T, const EVector& speciesNv) const
{
	return EVector::Zero(_numL);
}

EVector H2FromFiles::sinkv(double T, double n, const EVector& speciesNv) const
{
	return EVector::Zero(_numL);
}
