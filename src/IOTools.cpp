#include "IOTools.h"
#include "Error.h"
#include "global.h"

using namespace std;

ifstream IOTools::ifstreamFile(const string& file)
{
	ifstream input(file);
	if (!input)
	{
		string message = "File " + file + "not found.";
		Error::runtime(message);
	}
	DEBUG("Opened file " << file << endl);
	return input;
}

istringstream IOTools::istringstreamNextLine(ifstream& ifs)
{
	string line;
	getline(ifs, line);
	return istringstream(line);
}
