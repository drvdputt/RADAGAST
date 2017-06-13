#include "IOTools.h"
#include "global.h"

using namespace std;

ifstream IOTools::ifstreamFile(const string& file)
{
	ifstream input(file);
	if (!input)
	{
		string message = "File " + file + "not found.";
		cerr << message << endl;
		throw runtime_error(message);
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

