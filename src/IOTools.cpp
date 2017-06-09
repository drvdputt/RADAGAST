#include "IOTools.h"

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
	return input;
}

istringstream IOTools::istringstreamNextLine(ifstream& ifs)
{
	string line;
	getline(ifs, line);
	return istringstream(line);
}

