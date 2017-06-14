#include "HydrogenLevels.h"
#include "Constants.h"
#include "IOTools.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <iostream>
#include <vector>

using namespace std;

#define NLV 6

namespace
{
// commonly used indices
const int index2p = 1;
const int index2s = 2;

inline std::vector<int> twoJplus1range(int l)
{
	return l > 0 ? std::vector<int>({2 * l, 2 * l + 2}) : std::vector<int>({2});
}
}

HydrogenLevels::HydrogenLevels() : NLevel(makeNLv(), makeEv(), makeGv(), makeAvv(), makeExtraAvv())
{
	DEBUG("Constructed HydrogenLevels (without frequency grid)" << endl);
#ifndef HYDROGENLEVELS_HARDCODE
	readData();
#endif
}

HydrogenLevels::HydrogenLevels(const Array& frequencyv)
                : NLevel(frequencyv, makeNLv(), makeEv(), makeGv(), makeAvv(), makeExtraAvv())
{
	/* The correct way of grouping the l-levels is described in Hazy II, 6.11. When grouping the
	 initial levels of a collection of transitions, the collision strengths (bigUpsilon) need to
	 be summed, and transition probabilities (einstein Aij) need to be averaged over. Combining
	 the coefficients for different final levels comes down to
	 just summing them. */
	DEBUG("Constructed HydrogenLevels" << endl);
}

int HydrogenLevels::makeNLv() const
{
#ifdef HYDROGENLEVELS_HARDCODE
	return NLV;
#else
	return _chiantiNumLvl;
#endif
}

Eigen::VectorXd HydrogenLevels::makeEv() const
{
#ifdef HYDROGENLEVELS_HARDCODE
	Eigen::VectorXd the_ev(NLV);
	the_ev << 0., 82258.9191133, 82258.9543992821, 97492.304, 102823.904, 105291.657;
	// Energy in cm^-1, so multiply with hc
	return the_ev * Constant::PLANCKLIGHT;
#else

#endif
}

Eigen::VectorXd HydrogenLevels::makeGv() const
{
	Eigen::VectorXd the_gv(NLV);
	the_gv << 1, 1, 3, 9, 16, 25;
	return the_gv;
}

Eigen::MatrixXd HydrogenLevels::makeAvv() const
{
	/* Using NIST data
	 * (http://physics.nist.gov/cgi-bin/ASD/lines1.pl?unit=1&line_out=0&bibrefs=1&show_obs_wl=1&show_calc_wl=1&A_out=0&intens_out=1&allowed_out=1&forbid_out=1&conf_out=1&term_out=1&enrg_out=1&J_out=1&g_out=0&spectra=H%20I)
	 */

	/* Explicit implementation, for instructive purposes */

	// Decays to 1
	double A2p1 = 6.2649e8;
	double A2s1 = 2.495e-6;
	double A31 = 5.5751e7;
	double A41 = 1.2785e+07;
	double A51 = 4.1250e+06;

	// Decays to 2
	// l has to change by +-1, so 3 --> 2s can only happen from 3p
	// Two different initial states --> average
	double A3p2s = (2 * 2.2449e+07 + 4 * 2.2448e+07) / 6.;
	double A32s = (1 * 0. + 3 * A3p2s + 5 * 0.) / 9.;

	// Two different final states --> sum
	double A3s2p = 4.2097e+06 + 2.1046e+06;
	// Different final states and initial states --> both sum and average
	double A3d2p = (4 * 1.0775e+07 + 6 * 6.4651e+07 + 4 * 5.3877e+07) / 10.;
	double A32p = (1 * A3s2p + 3 * 0. + 5 * A3d2p) / 9.;

	// Keep this as a reference. The sum of A32s and A32p should equal this.
	double A32 = 4.4101e+07;
	printf("A32 = %e, A32p + A32s = %e\n", A32, A32p + A32s);

	double A4p2s = (2 * 9.6683e+06 + 4 * 9.6680e+06) / 6.;
	double A42s = 3 * A4p2s / 16.;

	double A4s2p = 1.7190e+06 + 8.5941e+05;
	double A4d2p = (4 * 3.4375e+06 + 6 * 2.0625e+07 + 4 * 1.7188e+07) / 10.;
	double A42p = (1 * A4s2p + 3 * 0. + 5 * A4d2p + 7 * 0.) / 16.;

	double A42 = 8.4193e+06;
	printf("A42 = %e, A42p + A42s = %e\n", A42, A42p + A42s);

	double A5p2s = (2 * 4.9484e+06 + 4 * 4.9483e+06) / 6.;
	double A52s = 3 * A5p2s / 25.;

	double A5s2p = 4.2955e+05 + 8.5920e+05;
	double A5d2p = (4 * 1.5709e+06 + 6 * 9.4254e+06 + 4 * 7.8548e+06) / 10.;
	double A52p = (A5s2p + 5 * A5d2p) / 25.;

	double A52 = 2.5304e+06;
	printf("A52 = %e, A52p + A52s = %e\n", A52, A52p + A52s);

	// Decays to 3
	double A43 = 8.9860e+06;
	double A53 = 2.2008e+06;

	// Decays to 4
	double A54 = 2.6993e+06;

	Eigen::MatrixXd the_avv(NLV, NLV);
	// clang-format off
	the_avv << 0, 0, 0, 0, 0, 0,
	           A2p1, 0, 0, 0, 0, 0,
	           A2s1, 0, 0, 0, 0, 0,
	           A31, A32p, A32s, 0, 0, 0,
	           A41, A42p, A42s, A43, 0, 0,
	           A51, A52p, A52s, A53, A54, 0;
	// clang-format on
	return the_avv;
}

Eigen::MatrixXd HydrogenLevels::makeExtraAvv() const
{
	/* Two photon continuum adds extra spontaneous decay rate to A2s1, but does not produce line
	 * radiation */
	Eigen::MatrixXd the_extraAvv(NLV, NLV);
	// clang-format off
	the_extraAvv << 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			8.23, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0;
	// clang-format on
	return the_extraAvv;
}

Eigen::MatrixXd HydrogenLevels::prepareCollisionMatrix(double T, double electronDensity,
                                                       double protonDensity) const
{
	Eigen::MatrixXd Cvv = Eigen::MatrixXd::Zero(NLV, NLV);

	auto fillInElectronCollisionRate = [&](size_t upper, size_t lower, double bigUpsilon) {
		double kT = Constant::BOLTZMAN * T;
		// Equation 6.17 of Hazy II (6.6 Collision strengths)
		Cvv(upper, lower) = bigUpsilon * 8.6291e-6 / gv(upper) / sqrt(T) * electronDensity;
		Cvv(lower, upper) = Cvv(upper, lower) * gv(upper) / gv(lower) *
		                    exp((ev(lower) - ev(upper)) / kT);
	};

	// data from 2002-anderson
	// Electron temperatures in electron volt
	const vector<double> grid_eT_eVv = {.5, 1., 3., 5., 10., 15., 20., 25.};

	double eT_eV = Constant::BOLTZMAN * T * Constant::ERG_EV;

	// Naively inter- and extrapolate this data linearly
	size_t iRight = TemplatedUtils::index(eT_eV, grid_eT_eVv);
	if (iRight == 0)
		iRight = 1;
	else if (iRight == grid_eT_eVv.size())
		iRight -= 1;
	size_t iLeft = iRight - 1;

	// Effective collision strength
	vector<double> bigUpsilon2s1v = {2.6e-1,  2.96e-1, 3.26e-1, 3.39e-1,
	                                 3.73e-1, 4.06e-1, 4.36e-1, 4.61e-1};
	double bigUpsilon2s1 = TemplatedUtils::interpolateLinear(
	                eT_eV, grid_eT_eVv[iLeft], grid_eT_eVv[iRight], bigUpsilon2s1v[iLeft],
	                bigUpsilon2s1v[iRight]);
	bigUpsilon2s1 = max(bigUpsilon2s1, 0.);
	fillInElectronCollisionRate(2, 0, bigUpsilon2s1);

	vector<double> bigUpsilon2p1v = {4.29e-1, 5.29e-01, 8.53e-01, 1.15e00,
	                                 1.81e00, 2.35e00,  2.81e00,  3.20e00};
	double bigUpsilon2p1 = TemplatedUtils::interpolateLinear(
	                eT_eV, grid_eT_eVv[iLeft], grid_eT_eVv[iRight], bigUpsilon2p1v[iLeft],
	                bigUpsilon2p1v[iRight]);
	bigUpsilon2p1 = max(bigUpsilon2s1, 0.);
	fillInElectronCollisionRate(1, 0, bigUpsilon2p1);

	// Important for two-photon continuum vs lyman is the l-changing collisions between 2s and
	// 2p 1964-Pengelly eq 43, assuming for qnl = qnl->nl' if l can only be l=1 or l=0
	double mu_m = Constant::HMASS_CGS / 2 / Constant::ELECTRONMASS;
	double constfactor = 9.93e-6 * sqrt(mu_m);

	//(6n^2(n^2 - l^2 - l - 1)) (eq 44)
	double D2p = 24;
	// even though the second term is zero
	double A2p = avv().row(index2p).sum() + extraAvv().row(index2p).sum();
	// (45 should be the correct one for DeltaE <<<)
	double twolog10Rc = 10.95 + log10(T / A2p / A2p / mu_m);
	double q2p2s = constfactor * D2p / sqrt(T) * (11.54 + log10(T / D2p / mu_m) + twolog10Rc);
	Cvv(index2p, index2s) = protonDensity * q2p2s;

	double D2s = 24 * 3;
	double A2s = avv().row(index2s).sum() + extraAvv().row(index2s).sum();
	twolog10Rc = 10.95 + log10(T / A2s / A2s / mu_m);
	double q2s2p = constfactor * D2s / sqrt(T) * (11.54 + log10(T / D2s / mu_m) + twolog10Rc);
	Cvv(index2s, index2p) = protonDensity * q2s2p;

	return Cvv;
}

Array HydrogenLevels::boundBoundContinuum(const Solution& s) const
{
	Array result(frequencyv().size());
	// 1984-Nussbaumer
	double constFactor = Constant::PLANCK / Constant::FPI * s.nv(index2s);
	double nu0 = (ev(index2s) - ev(0)) / Constant::PLANCK;
	double C = 202.0; // s-1
	double alpha = .88;
	double beta = 1.53;
	double gam = .8;
	for (size_t iFreq = 0; frequencyv()[iFreq] < nu0; iFreq++)
	{
		double y = frequencyv()[iFreq] / nu0;
		double y1miny = y * (1 - y);
		double pow4y1miny_gam = pow(4 * y1miny, gam);
		double Py = C * (y1miny * (1 - pow4y1miny_gam) +
		                 alpha * pow(y1miny, beta) * pow4y1miny_gam);
		result[iFreq] = constFactor * y * Py;
	}
	return result;
}

void HydrogenLevels::readData()
{
	if (_dataReady)
		return;
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
	_chiantiAvv.resize(_chiantiNumLvl, _chiantiNumLvl);

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
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

		// Zero means two-photon transition, see CHIANTI user guide
		_chiantiAvv(upperIndex - 1, lowerIndex) = wavAngstrom > 0 ? A : 0;

		getline(einsteinFile, line);
	}
	einsteinFile.close();

	//---------------------//
	// READ COLLISION DATA //
	//---------------------//
	int andersonIndex = 1;
	for (int n = 0; n <= 5; n++)
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
	_dataReady = true;
	exit(1);
}

double HydrogenLevels::energy(int n, int l) const
{
	// Take an average over the j states
	double esum = 0;
	for (int twoJplus1 : twoJplus1range(l))
		esum += _chiantiLevelv[indexCHIANTI(n, l, twoJplus1)].e * twoJplus1;
	return esum / (4 * l + 2);
}

double HydrogenLevels::energy(int n) const
{
	// Average over the l states
	double esum = 0;
	for (int l = 0; l < n; l++)
		esum += energy(n, l) * (2 * l + 1);
	return esum / n*n;
}

double HydrogenLevels::einsteinA(int ni, int li, int nf, int lf) const
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

double HydrogenLevels::einsteinA(int ni, int li, int nf) const
{
	// sum over the final l states
	double Asum = 0;
	for (int lf = 0; lf < nf; lf++)
		Asum += einsteinA(ni, li, nf, lf);
	return Asum;
}

double HydrogenLevels::einsteinA(int ni, int nf) const
{
	// average over the initial l states
	double Asum = 0;
	for (int li = 0; li < ni; li++)
	{
		Asum += einsteinA(ni, li, nf) * (2 * li + 1);
	}
	// divide by multiplicity of ni state to get the average (factor two has been dropped in
	// enumerator and denominator)
	return Asum / (ni * ni);
}

double HydrogenLevels::eCollisionStrength(int ni, int li, int nf, int lf, double T_eV) const
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

double HydrogenLevels::eCollisionStrength(int ni, int nf, int lf, double T_eV) const
{
	// Collision strengths must be summed over all initial states
	double Upsilonsum = 0;
	for (int li = 0; li < ni; li++)
		Upsilonsum += eCollisionStrength(ni, li, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenLevels::eCollisionStrength(int ni, int nf, double T_eV) const
{
	// Also sum over all final states
	double Upsilonsum = 0;
	for (int lf = 0; lf < nf; lf++)
		Upsilonsum += eCollisionStrength(ni, nf, lf, T_eV);
	return Upsilonsum;
}
