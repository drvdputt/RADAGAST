#include "IOTools.h"
#include "Error.h"
#include "global.h"

using namespace std;

ifstream IOTools::ifstreamFile(const string& file)
{
	ifstream input(file);
	if (!input)
		Error::runtime("Input file " + file + "not found.");
	DEBUG("Opened file (read)" << file << endl);
	return input;
}

ofstream IOTools::ofstreamFile(const string& file)
{
	ofstream output(file);
	if (!output.is_open())
		Error::runtime("Output file " + file + " could not be opened.");
	DEBUG("Opened file (write)" << file << endl);
	return output;
}

istringstream IOTools::istringstreamNextLine(ifstream& ifs)
{
	string line;
	getline(ifs, line);
	return istringstream(line);
}
