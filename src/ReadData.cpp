#include "Constants.h"
#include "ReadData.h"

#include <fstream>
#include <iterator>
#include <sstream>

using namespace std;

void ReadData::recombinationContinuum(string file, vector<double>& fileFrequencyv, vector<double>& fileThresholdv,
		vector<double>& fileTemperaturev, vector<vector<double>>& fileGammaDaggervv)
{
	ifstream input(file);
	if (!input)
		throw std::runtime_error("File " + file + " not found.");

	size_t numcol, numrow;

	// The line number will not count any lines starting with #
	int lineNr = 0;
	string line;
	while (getline(input, line))
	{
		istringstream iss(line);
		if (iss.peek() != '#')
		{
			if (lineNr == 0)
			{
				iss >> numcol >> numrow;
				fileGammaDaggervv.resize(numrow, std::vector<double>(numcol, 0.));
			}
			if (lineNr == 1)
			{
				double log10T;
				while (iss >> log10T)
					fileTemperaturev.push_back(log10T);
			}
			if (lineNr > 1)
			{
				int flag;
				double energy;
				iss >> flag >> energy;

				double frequency = energy * Constant::RYDBERG / Constant::PLANCK;
				fileFrequencyv.push_back(frequency);
				if (flag)
					fileThresholdv.push_back(frequency);

				auto it = fileGammaDaggervv[lineNr - 2].begin();
				double gammaDagger;
				while (iss >> gammaDagger)
				{
					*it = gammaDagger;
					it++;
				}
			}
			lineNr++;
		}
	}
}
