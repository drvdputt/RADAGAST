#include "HydrogenLevels.h"
#include "Constants.h"
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
}

HydrogenLevels::HydrogenLevels() : NLevel(makeNLv(), makeEv(), makeGv(), makeAvv(), makeExtraAvv())
{
	DEBUG("Constructed HydrogenLevels (without frequency grid)" << endl);
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

int HydrogenLevels::makeNLv() const { return NLV; }

Eigen::VectorXd HydrogenLevels::makeEv() const
{
	Eigen::VectorXd the_ev(NLV);
	the_ev << 0., 82258.9191133, 82258.9543992821, 97492.304, 102823.904, 105291.657;
	// Energy in cm^-1, so multiply with hc
	return the_ev * Constant::PLANCKLIGHT;
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
