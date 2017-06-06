#include "TwoLevel.h"
#include "Constants.h"
#include "Sanity.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <algorithm>
#include <exception>
#include <vector>

TwoLevel::TwoLevel(const Array& frequencyv)
                : NLevel(frequencyv, makeNLv(), makeEv(), makeGv(), makeAvv(), makeExtraAvv())
{
}

int TwoLevel::makeNLv() const
{
	/* toy model of CII 158 um from https://www.astro.umd.edu/~jph/N-level.pdf bottom of
	 page 4. This makes collisional deexcitation important for densities > 20-50 / cm3. */
	// Level and transition information
	//	_Ev << 0, .00786 / Constant::ERG_EV;
	//	_gv << 1, 1;
	//	_nv << _n, 0;

	// The only spontaneous transition is A_10
	//	_Avv << 0, 0,
	//		   2.29e-6, 0;

	return 2;
}

Eigen::VectorXd TwoLevel::makeEv() const
{
	Eigen::VectorXd the_ev(2);
	the_ev << 0, .00786 / Constant::ERG_EV;
	return the_ev;
}

Eigen::VectorXd TwoLevel::makeGv() const
{
	Eigen::VectorXd the_gv(2);
	the_gv << 1, 1;
	return the_gv;
}

Eigen::MatrixXd TwoLevel::makeAvv() const
{
	Eigen::MatrixXd the_avv(2, 2);
	// clang-format off
	the_avv << 0, 0,
		   2.29e-6, 0;
	// clang-format on
	return the_avv;
}

Eigen::MatrixXd TwoLevel::prepareCollisionMatrix(double T, double electronDensity, double protonDensity) const
{
	Eigen::MatrixXd Cvv = Eigen::MatrixXd::Zero(2, 2);

	// Need separate contributions for number of protons and electrons
	// Toy implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
	// is actually for electron collisions only, but let's treat all collision partners
	// this way for now
	double beta = 8.629e-6;

	// also take some values from the bottom of page 4
	// Gamma = 2.15 at 10000 K and 1.58 at 1000 K
	double bigUpsilon10 = (T - 1000) / 9000 * 2.15 + (10000 - T) / 9000 * 1.58;

	Cvv(1, 0) = beta / sqrt(T) * bigUpsilon10 / gv(1);
	Cvv(0, 1) = Cvv(1, 0) * gv(1) / gv(0) * exp(-(ev(1) - ev(0)) / Constant::BOLTZMAN / T);
	return Cvv;
}
