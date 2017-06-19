#include "HydrogenFromFiles.h"

#include "Constants.h"
#include "Error.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "global.h"

#include <cassert>
#include <fstream>
#include <sstream>

using namespace std;

namespace
{
const int nMax = 5;

inline std::vector<int> twoJplus1range(int l)
{
	// If l > 0, then j can be either l + 0.5 or l - 0.5. For an s state, j is always 1/2 and
	// 2j+1 can only be 2
	return l > 0 ? std::vector<int>({2 * l, 2 * l + 2}) : std::vector<int>({2});
}
}

HydrogenFromFiles::HydrogenFromFiles(int resolvedUpTo) : _resolvedUpTo(resolvedUpTo)
{
	if (_resolvedUpTo > nMax)
		Error::runtime("There are currently only five shells for the hydrogen model");
	// Read-in, and steps that are safe during read-in
	readData();
	// Steps that need to happen after read-in
	prepareForOutput();
}

void HydrogenFromFiles::readData()
{
	const std::string basename(REPOROOT "/dat/CHIANTI_8.0.6_data/h/h_1/h_1");

	//-----------------//
	// READ LEVEL DATA //
	//-----------------//
	ifstream elvlc = IOTools::ifstreamFile(basename + ".elvlc");
	string line;
	getline(elvlc, line);
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

		int n;
		istringstream(config) >> n;    // Get the first number from the config string,
		int l = _lNumberm.at(lSymbol); // Translate the angular momentum letter
		int twoJplus1 = static_cast<int>(2 * j + 1); // Store 2j+1 (static cast to make it
		                                             // clear that j is not integer)
		double e = observedEnergy * Constant::LIGHT *
		           Constant::PLANCK; // Convert from cm-1 to erg
		DEBUG("n " << n << " l " << l << " 2j+1 " << twoJplus1 << " e "
		           << e * Constant::ERG_EV << endl);

		// The level indices in the data structures will go from 0 to number of levels minus
		// one
		// The quantum numbers are also used as keys in a map, so we can quickly retrieve
		// the index for a given configuration.
		// The level indices in the file and the map will go from 1 to the number of levels
		_chiantiLevelv.emplace_back(n, l, twoJplus1, e);
		_nljToChiantiIndexm.insert({{n, l, twoJplus1}, lvIndex - 1});

		getline(elvlc, line);
	}
	elvlc.close();
	_chiantiNumLvl = _chiantiLevelv.size();

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
	_chiantiAvv = Eigen::MatrixXd::Zero(_chiantiNumLvl, _chiantiNumLvl);
	ifstream wgfa = IOTools::ifstreamFile(basename + ".wgfa");
	getline(wgfa, line);
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

		getline(wgfa, line);
	}
	wgfa.close();

	DEBUG("All chianti A coefficients:" << endl);
	DEBUG(_chiantiAvv << endl);

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

	ifstream h_coll_str = IOTools::ifstreamFile(REPOROOT "/dat/h_coll_str.dat");
	getline(h_coll_str, line);
	getline(h_coll_str, line);
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
		getline(h_coll_str, line);
	}
	h_coll_str.close();
}

void HydrogenFromFiles::prepareForOutput()
{
	//-----------------------//
	// SET OUTPUT PARAMETERS //
	//-----------------------//
	_numL = 0;
	_levelOrdering.clear();
	int n = 1;
	while (n <= _resolvedUpTo)
	{
		for (int l = 0; l < n; l++)
		{
			_numL++;
			_levelOrdering.emplace_back(n, l, energy(n, l));
		}
		n++;
	}
	while (n <= nMax)
	{
		_numL++;
		_levelOrdering.emplace_back(n, energy(n));
		n++;
	}
}

int HydrogenFromFiles::numLv() const { return _numL; }

Eigen::VectorXd HydrogenFromFiles::ev() const
{
	Eigen::VectorXd the_ev(_numL);
	for (int i = 0; i < _numL; i++)
	{
		const LevelInfo& lvInfo = _levelOrdering[i];
		// l = -1 means collapsed, so:
		the_ev[i] = lvInfo.e();
	}
	DEBUG("energy vector" << endl);
	DEBUG(the_ev << endl);
	return the_ev;
}

Eigen::VectorXd HydrogenFromFiles::gv() const
{
	Eigen::VectorXd the_gv(_numL);
	for (int i = 0; i < _numL; i++)
	{
		const LevelInfo& lvInfo = _levelOrdering[i];
		// collapsed ? n^2 : 2l+1
		the_gv[i] = lvInfo.l() < 0 ? lvInfo.n() * lvInfo.n() : 2 * lvInfo.l() + 1;
	}
	DEBUG("degeneracy vector" << endl);
	DEBUG(the_gv << endl);
	return the_gv;
}

Eigen::MatrixXd HydrogenFromFiles::avv() const
{
	Eigen::MatrixXd the_avv(_numL, _numL);
	for (int i = 0; i < _numL; i++)
	{
		const LevelInfo& initial = _levelOrdering[i];
		int ni = initial.n();
		int li = initial.l();
		for (int f = 0; f < _numL; f++)
		{
			const LevelInfo& final = _levelOrdering[f];
			int nf = final.n();
			int lf = final.l();
			// Both resolved
			if (li >= 0 && lf >= 0)
				the_avv(i, f) = einsteinA(ni, li, nf, lf);
			// Collapsed to resolved
			else if (li < 0 && lf >= 0)
				the_avv(i, f) = einsteinA(ni, nf, lf);
			// Resolved to collapsed. The collapsed-collapsed equivalent should always
			// be 0 in this case.
			else if (li >= 0 && lf < 0)
			{
				the_avv(i, f) = einsteinA(ni, nf);
				assert(the_avv(i, f) == 0.);
			}
			// Collapsed to collapsed
			else
				the_avv(i, f) = einsteinA(ni, nf);
		}
	}
	DEBUG("Einstein A matrix:" << endl);
	DEBUG(the_avv << endl);
	return the_avv;
}

Eigen::MatrixXd HydrogenFromFiles::extraAvv() const
{
	Eigen::MatrixXd the_extra = Eigen::MatrixXd::Zero(_numL, _numL);
	int index1s = -1;
	int index2p = -1;
	int i = 0;
	while (index1s < 0 || index2p < 0)
	{
		const LevelInfo& lvInfo = _levelOrdering[i];
		if (lvInfo.n() == 1)
			index1s = i;
		if (lvInfo.n() == 2 && lvInfo.l() == 1)
			index2p = i;
		if (i >= _numL)
			Error::runtime("Can't find 2p level for hydrogen");
		i++;
	}
	// Hardcoded the two-photon decay of 2p to 1s
	the_extra(index2p, index1s) = 8.229;
	DEBUG("Extra A" << endl);
	DEBUG(the_extra << endl);
	return the_extra;
}

Eigen::MatrixXd HydrogenFromFiles::cvv(double T, double ne, double np) const
{
	// Calculate the temperature in erg and in electron volt
	double kT = Constant::BOLTZMAN * T;
	double T_eV = kT * Constant::ERG_EV;

	Eigen::MatrixXd the_cvv = Eigen::MatrixXd::Zero(_numL, _numL);
	// Electron contributions (n-changing)
	for (int i = 0; i < _numL; i++)
	{
		const LevelInfo& ini = _levelOrdering[i];
		for (int f = 0; f < _numL; f++)
		{
			const LevelInfo& fin = _levelOrdering[f];
			// for downward transitions, calculate the collision rate, and derive the
			// rate for the corresponding upward transition too
			if (ini.e() > fin.e())
			{
				double UpsilonDown = eCollisionStrength(ini, fin, T_eV);
				double Cif = UpsilonDown * 8.6291e-6 / ini.g() / sqrt(T) * ne;
				double Cfi = Cif * ini.g() / fin.g() * exp(fin.e() - ini.e() / kT);
				the_cvv(i, f) += Cif;
				the_cvv(f, i) += Cfi;
			}
			// for upward transitions do nothing because we already covered them above
		}
	}

	// TODO: Proton contributions (l-changing)

	return the_cvv;
}

double HydrogenFromFiles::energy(int n, int l) const
{
	// Take an average over the j states
	double esum = 0;
	for (int twoJplus1 : twoJplus1range(l))
		esum += _chiantiLevelv[indexCHIANTI(n, l, twoJplus1)].e() * twoJplus1;
	return esum / (4 * l + 2);
}

double HydrogenFromFiles::energy(int n) const
{
	// Average over the l states
	double esum = 0;
	for (int l = 0; l < n; l++)
		esum += energy(n, l) * (2 * l + 1);
	return esum / n * n;
}

double HydrogenFromFiles::einsteinA(int ni, int li, int nf, int lf) const
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

double HydrogenFromFiles::einsteinA(int ni, int nf, int lf) const
{
	// average over the initial l states
	double Asum = 0;
	for (int li = 0; li < ni; li++)
		Asum += einsteinA(ni, li, nf, lf) * (2 * li + 1);
	return Asum / (ni * ni);
}

double HydrogenFromFiles::einsteinA(int ni, int nf) const
{
	// sum over the final l states
	double Asum = 0;
	for (int lf = 0; lf < nf; lf++)
		Asum += einsteinA(ni, nf, lf);
	return Asum;
}

double HydrogenFromFiles::eCollisionStrength(int ni, int li, int nf, int lf, double T_eV) const
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
		Error::runtime("This function should only be used for downward collisional "
		               "transitions");

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

double HydrogenFromFiles::eCollisionStrength(int ni, int nf, int lf, double T_eV) const
{
	// Collision strengths must be summed over all initial states
	double Upsilonsum = 0;
	for (int li = 0; li < ni; li++)
		Upsilonsum += eCollisionStrength(ni, li, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenFromFiles::eCollisionStrength(int ni, int nf, double T_eV) const
{
	// Also sum over all final states
	double Upsilonsum = 0;
	for (int lf = 0; lf < nf; lf++)
		Upsilonsum += eCollisionStrength(ni, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenFromFiles::eCollisionStrength(const LevelInfo& initial, const LevelInfo& final,
                                             double T_eV) const
{
	// No output for upward transitions
	if (initial.e() < final.e())
		return 0;
	else
	{
		// Resolved-resolved
		if (initial.l() >= 0 && final.l() >= 0)
			return eCollisionStrength(initial.n(), initial.l(), final.n(), final.l(),
			                          T_eV);
		// Collapsed-resolved
		else if (initial.l() < 0 && final.l() >= 0)
			return eCollisionStrength(initial.n(), final.n(), final.l(), T_eV);
		// Resolved-collapsed (should not be called, and if it is, the collapsed-collapsed
		// result should be zero)
		else if (initial.l() >= 0 && final.l() < 0)
		{
			double eCollStr = eCollisionStrength(initial.n(), final.n(), T_eV);
			assert(eCollStr == 0);
			return eCollStr;
		}
		// Collapsed-collapsed
		else
			return eCollisionStrength(initial.n(), final.n(), T_eV);
	}
}
