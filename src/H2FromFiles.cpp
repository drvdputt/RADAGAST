#include "H2FromFiles.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"

// TODO: make model

using namespace std;

H2FromFiles::H2FromFiles() { readData(); }

size_t H2FromFiles::numLv() const { return 1; }

size_t H2FromFiles::indexOutput(ElectronicState eState, int j, int v) const
{
	// This cast is safe, as the default underlying type of an enum is int.
	return _evjToIndex.at({static_cast<int>(eState), j, v});
}

void H2FromFiles::readData()
{
	const string location{REPOROOT "/dat/h2"};

	//-----------------//
	// READ LEVEL DATA //
	//-----------------//
	ifstream energyX = IOTools::ifstreamFile(location + "/energy_X.dat");
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
	DEBUG("Read in " << _levelv.size() << " H2 levels" << endl);
	energyX.close();

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
	_avv = EMatrix::Zero(_levelv.size(), _levelv.size());

	ifstream transprobX = IOTools::ifstreamFile(location + "/transprob_X.dat");
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
#endif
}

EVector H2FromFiles::ev() const { return EVector::Zero(1); }

EVector H2FromFiles::gv() const { return EVector::Constant(1, 1); }

EMatrix H2FromFiles::avv() const { return EMatrix::Zero(1, 1); }

EMatrix H2FromFiles::extraAvv() const { return EMatrix::Zero(1, 1); }

EMatrix H2FromFiles::cvv(double T, double ne, double np) const { return EMatrix::Zero(1, 1); }

EVector H2FromFiles::sourcev(double T, double ne, double np) const { return EVector::Zero(1); }

EVector H2FromFiles::sinkv(double T, double n, double ne, double np) const
{
	return EVector::Zero(1);
}
