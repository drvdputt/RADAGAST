#include "H2FromFiles.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"

// TODO: make model

using namespace std;

H2FromFiles::H2FromFiles() { readData(); }

int H2FromFiles::numLv() const { return 1; }

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
			/* Include only the levels with the "data missing" comment, not those that
			   are "#extra data". */
			// if (comment != "#data missing")
			// continue;
			DEBUG(word1 + word2);
		}

		// Get the numbers
		int v, j;
		iss >> v >> j;
		double e;
		iss >> e;

		// Create a new level object
		_levelv.emplace_back(H2FromFiles::ElectronicState::X, j, v, e);
	}
	DEBUG("Read in " << _levelv.size() << " H2 levels" << endl);
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
