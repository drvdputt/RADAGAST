#include "HydrogenDataProvider.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "global.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace
{
const int nMax = 5;

inline std::vector<int> twoJplus1range(int l)
{
	return l > 0 ? std::vector<int>({2 * l, 2 * l + 2}) : std::vector<int>({2});
}
}

HydrogenDataProvider::HydrogenDataProvider(int resolvedUpTo) : _resolvedUpTo(resolvedUpTo)
{
	if (_resolvedUpTo > nMax)
		throw "There are currently only five shells for the hydrogen model";
	readData();
}

int HydrogenDataProvider::numLv() const { return _numL; }

Eigen::VectorXd HydrogenDataProvider::ev() const
{
	Eigen::VectorXd the_ev(_numL);
	for (int i = 0; i < _numL; i++)
	{
		NLPair nl = _levelOrdering[i];
		// l = -1 means collapsed, so:
		the_ev[i] = nl.l < 0 ? energy(nl.n) : energy(nl.n, nl.l);
	}
	return the_ev;
}

Eigen::VectorXd HydrogenDataProvider::gv() const
{
	Eigen::VectorXd the_gv(_numL);
	for (int i = 0; i < _numL; i++)
	{
		NLPair nl = _levelOrdering[i];
		// collapsed ? n^2 : 2l+1
		the_gv[i] = nl.l < 0 ? nl.n * nl.n : 2 * nl.l + 1;
	}
	return the_gv;
}

Eigen::MatrixXd HydrogenDataProvider::avv() const
{
	Eigen::MatrixXd the_avv(_numL, _numL);
	for (int i = 0; i < _numL; i++)
	{
		NLPair nli = _levelOrdering[i];
		int ni = nli.n;
		int li = nli.l;
		for (int f = 0; f < _numL; f++)
		{
			NLPair nlf = _levelOrdering[f];
			int nf = nlf.n;
			int lf = nlf.l;
			// Both resolved
			if (li > 0 && lf > 0)
				the_avv(i, f) = einsteinA(ni, li, nf, lf);
			// Collapsed to resolved
			else if (li < 0 && lf > 0)
				the_avv(i, f) = einsteinA(ni, nf, lf);
			// Resolved to collapsed. The collapsed-collapsed equivalent should always
			// be 0 in this case.
			else if (li > 0 && lf < 0)
			{
				the_avv(i, f) = einsteinA(ni, nf);
				assert(the_avv(i, f) == 0.);
			}
			// Collapsed to collapsed
			else
				the_avv(i, f) = einsteinA(ni, nf);
		}
	}
	return the_avv;
}

Eigen::MatrixXd HydrogenDataProvider::extraAvv() const
{
	Eigen::MatrixXd the_extra = Eigen::MatrixXd::Zero(_numL, _numL);
	int index1s = -1;
	int index2p = -1;
	int i = 0;
	while (index1s < 0 || index2p < 0)
	{
		NLPair nl = _levelOrdering[i];
		if (nl.n == 1) index1s = i;
		if (nl.n == 2 && nl.l == 1) index2p = i;
		if (i >= _numL) throw "Can't find 2p level for hydrogen";
		i++;
	}
	// Hardcoded the two-photon decay of 2p to 1s
	the_extra(index1s, index2p) = 8.229;
	return the_extra;
}

void HydrogenDataProvider::readData()
{
	const std::string basename(REPOROOT "/dat/CHIANTI_8.0.6_data/h/h_1/h_1");

	//-----------------//
	// READ LEVEL DATA //
	//-----------------//
	ifstream levelFile = IOTools::ifstreamFile(basename + ".elvlc");
	string line;
	getline(levelFile, line);
	while (line.compare(1, 2, "-1"))
	{
		// Read the different parts of the line
		int lvIndex, twoSplus1;
		string config;
		char lSymbol;
		double j, observedEnergy, theoreticalEnergy;
		istringstream(line) >> lvIndex >> config >> twoSplus1 >> lSymbol >> j >>
		                observedEnergy >> theoreticalEnergy;
		DEBUG(lvIndex << " " << config << " " << twoSplus1 << " " << lSymbol << " " << j
		              << " " << observedEnergy << " " << theoreticalEnergy << " " << endl);

		levelInfo lvInfo;
		istringstream(config) >> lvInfo.n; // Get the first number from the config string,
		lvInfo.l = _lNumberm.at(lSymbol);  // Translate the angular momentum letter
		lvInfo.twoJplus1 = static_cast<int>(2 * j + 1); // Store 2j+1
		lvInfo.e = observedEnergy;
		DEBUG("n " << lvInfo.n << " l " << lvInfo.l << " 2j+1 " << lvInfo.twoJplus1 << " e "
		           << lvInfo.e << endl);

		// The level indices in the data structures will go from 0 to number of levels minus
		// one
		// The quantum numbers are also used as keys in a map, so we can quickly retrieve
		// the index for a given configuration.
		// The level indices in the file and the map will go from 1 to the number of levels
		_chiantiLevelv.push_back(lvInfo);
		_nljToChiantiIndexm.insert({{lvInfo.n, lvInfo.l, lvInfo.twoJplus1}, lvIndex - 1});

		getline(levelFile, line);
	}
	levelFile.close();
	_chiantiNumLvl = _chiantiLevelv.size();

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
	_chiantiAvv.resize(_chiantiNumLvl, _chiantiNumLvl);
	ifstream einsteinFile = IOTools::ifstreamFile(basename + ".wgfa");
	getline(einsteinFile, line);
	while (line.compare(1, 2, "-1"))
	{
		int leftIndex, rightIndex;
		double wavAngstrom, gf, A;
		istringstream(line) >> leftIndex >> rightIndex >> wavAngstrom >> gf >> A;

		// A comment in the cloudy code recommended to do this, as there are apparently some
		// files in the CHIANTI database where the left index represents the upper level of
		// the transition:
		int upperIndex = max(leftIndex, rightIndex);
		int lowerIndex = min(leftIndex, rightIndex);

		DEBUG(lowerIndex << " " << upperIndex << " " << wavAngstrom << " " << gf << " " << A
		                 << endl);

		// Zero means two-photon transition, see CHIANTI user guide.
		_chiantiAvv(upperIndex - 1, lowerIndex) = wavAngstrom > 0 ? A : 0;

		getline(einsteinFile, line);
	}
	einsteinFile.close();

	//---------------------//
	// READ COLLISION DATA //
	//---------------------//
	int andersonIndex = 1;
	for (int n = 0; n <= nMax; n++)
	{
		for (int l = 0; l < n; l++)
		{
			_nlToAndersonIndexm.insert({{n, l}, andersonIndex});
			andersonIndex++;
		}
	}

	ifstream andersonFile = IOTools::ifstreamFile(REPOROOT "/dat/h_coll_str.dat");
	getline(andersonFile, line);
	getline(andersonFile, line);
	while (line.compare(0, 2, "-1"))
	{
		istringstream iss = istringstream(line);
		int upperIndex, lowerIndex;
		iss >> upperIndex >> lowerIndex;
		Array Upsilonv(8);
		for (int t = 0; t < 8; t++)
		{
			DEBUG("temp" << _andersonTempv[t]);
			iss >> Upsilonv[t];
			DEBUG(" " << Upsilonv[t] << " ");
		}
		DEBUG(endl);
		_andersonUpsilonvm[{upperIndex, lowerIndex}] = Upsilonv;
		getline(andersonFile, line);
	}
	andersonFile.close();

	//--------------------------------//
	// SET THE TOTAL NUMBER OF LEVELS //
	//--------------------------------//
	_numL = 0;
	_levelOrdering.clear();
	int n = 1;
	while (n <= _resolvedUpTo)
	{
		for (int l = 0; l < n; l++)
		{
			_numL++;
			_levelOrdering.push_back({.n = n, .l = l});
		}
		n++;
	}
	while (n <= nMax)
	{
		_numL++;
		_levelOrdering.push_back({.n = n, .l = -1});
		n++;
	}
}

double HydrogenDataProvider::energy(int n, int l) const
{
	// Take an average over the j states
	double esum = 0;
	for (int twoJplus1 : twoJplus1range(l))
		esum += _chiantiLevelv[indexCHIANTI(n, l, twoJplus1)].e * twoJplus1;
	return esum / (4 * l + 2);
}

double HydrogenDataProvider::energy(int n) const
{
	// Average over the l states
	double esum = 0;
	for (int l = 0; l < n; l++)
		esum += energy(n, l) * (2 * l + 1);
	return esum / n * n;
}

double HydrogenDataProvider::einsteinA(int ni, int li, int nf, int lf) const
{
	if (ni < nf)
		return 0.;
	else
	{
		double Asum = 0;
		// sum over the final j states
		for (int twoJplus1f : twoJplus1range(lf))
		{
			int lIndex = indexCHIANTI(nf, lf, twoJplus1f);

			// average over the initial j states
			for (int twoJplus1i : twoJplus1range(li))
			{
				int uIndex = indexCHIANTI(ni, li, twoJplus1i);
				Asum += _chiantiAvv(uIndex, lIndex) * twoJplus1i;
			}
		}
		// divide by multiplicity of li state to get the average
		return Asum / (4 * li + 2);
	}
}

double HydrogenDataProvider::einsteinA(int ni, int nf, int lf) const
{
	// average over the initial l states
	double Asum = 0;
	for (int li = 0; li < ni; li++)
		Asum += einsteinA(ni, li, nf, lf) * (2 * li + 1);
	return Asum / (ni * ni);
}

double HydrogenDataProvider::einsteinA(int ni, int nf) const
{
	// sum over the final l states
	double Asum = 0;
	for (int lf = 0; lf < nf; lf++)
		Asum += einsteinA(ni, nf, lf);
	return Asum;
}

double HydrogenDataProvider::eCollisionStrength(int ni, int li, int nf, int lf, double T_eV) const
{
	// Alternatively to all the mappy things below, we could use a nested vector filled up with
	// zeros
	auto indexMapEnd = _nlToAndersonIndexm.end();
	// When the level is not included in the range of the data, the result is 0
	auto uIndexIt = _nlToAndersonIndexm.find({ni, li});
	auto lIndexIt = _nlToAndersonIndexm.find({nf, lf});
	if (uIndexIt == indexMapEnd || lIndexIt == indexMapEnd)
		return 0;

	int uIndex = uIndexIt->second;
	int lIndex = lIndexIt->second;
	if (uIndex <= lIndex)
		throw "This function should only be used for downward collisional transitions ";

	// When the levels are included, but the specific transition isn't, the result is also zero
	auto UpsilonvIt = _andersonUpsilonvm.find({uIndex, lIndex});
	if (UpsilonvIt == _andersonUpsilonvm.end())
		return 0;
	else
	{
		// If data is available for this transition, interpolate at the given temperature
		// (in eV)
		double Upsilon_ip = TemplatedUtils::evaluateLinInterpf<double, Array, Array>(
		                T_eV, _andersonTempv, UpsilonvIt->second);
		return Upsilon_ip;
	}
}

double HydrogenDataProvider::eCollisionStrength(int ni, int nf, int lf, double T_eV) const
{
	// Collision strengths must be summed over all initial states
	double Upsilonsum = 0;
	for (int li = 0; li < ni; li++)
		Upsilonsum += eCollisionStrength(ni, li, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenDataProvider::eCollisionStrength(int ni, int nf, double T_eV) const
{
	// Also sum over all final states
	double Upsilonsum = 0;
	for (int lf = 0; lf < nf; lf++)
		Upsilonsum += eCollisionStrength(ni, nf, lf, T_eV);
	return Upsilonsum;
}
