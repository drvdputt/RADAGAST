#include "RecombinationRate.hpp"
#include "IOTools.hpp"
#include "TemplatedUtils.hpp"

using namespace std;

const map<char, int> spdf_to_l = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3}, {'G', 4}, {'H', 5}, {'I', 6}, {'J',7}};

HydrogenADF48::HydrogenADF48()
{
	readADF48File("dat/adf48-cut");
}

HydrogenADF48::~HydrogenADF48() = default;

void HydrogenADF48::readADF48File(const string& path)
{
	auto ifs = IOTools::ifstreamRepoFile(path);
	string line;

	// indexing
	getline(ifs, line);
	int index = 0;
	while(getline(ifs, line))
	{
		if (line.empty())
			break;
		int i, g;
		string nl_string;
		istringstream(line) >> i >> g >> nl_string;
		int n = nl_string[0] - '0';
		int l = spdf_to_l.at(nl_string[1]);
		_nlToIndexm.insert({{n,l}, index});
		index++;
	}
	// temperature
	getline(ifs, line);
	getline(ifs, line);
	auto tv = IOTools::allNumbersFromNextLine(line);
	_temperaturev = Array(tv.data(), tv.size());

	// recombination coefficient
	getline(ifs,line);
	while(getline(ifs, line))
	{
		vector<double> alphav = IOTools::allNumbersFromNextLine(line);
		Error::equalCheck("Number of columns in adf48 file", alphav.size(), tv.size());
		_alphavv.emplace_back(alphav);
	}
}

double HydrogenADF48::alpha(int n, int l, double T) const
{
	// Grab the right line
	const vector<double>& alphav = _alphavv[_nlToIndexm.at({n, l})];

	// Interpolate on the temperature
	return TemplatedUtils::evaluateLinInterpf(T, _temperaturev, alphav);
}
