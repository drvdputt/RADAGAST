#include "TwoLevelHardcoded.h"
#include "Constants.h"

TwoLevelHardcoded::TwoLevelHardcoded()
{
	/* toy model of CII 158 um from https://www.astro.umd.edu/~jph/N-level.pdf bottom of page
	   4. This makes collisional deexcitation important for densities > 20-50 / cm3. */
	// Level and transition information
	//	_Ev << 0, .00786 / Constant::ERG_EV;
	//	_gv << 1, 1;
	//	_nv << _n, 0;

	// The only spontaneous transition is A_10
	//	_Avv << 0, 0,
	//		   2.29e-6, 0;
	the_ev << 0, .00786 / Constant::ERG_EV;
	the_gv << 1, 1;
}

int TwoLevelHardcoded::numLv() const { return 2; }

Eigen::VectorXd TwoLevelHardcoded::ev() const { return the_ev; }

Eigen::VectorXd TwoLevelHardcoded::gv() const { return the_gv; }

Eigen::MatrixXd TwoLevelHardcoded::avv() const
{
	Eigen::Matrix2d the_avv;
	// clang-format off
	the_avv << 0, 0,
		   2.29e-6, 0;
	// clang-format on
	return the_avv;
}

Eigen::MatrixXd TwoLevelHardcoded::extraAvv() const { return Eigen::Matrix2d::Zero(); }

Eigen::MatrixXd TwoLevelHardcoded::cvv(double T, double /* unused ne */,
                                       double /* unused np */) const
{
	Eigen::MatrixXd Cvv = Eigen::MatrixXd::Zero(2, 2);

	/* Need separate contributions for number of protons and electrons. Toy implementation below,
	   inspired by https://www.astro.umd.edu/~jph/N-level.pdf is actually for electron
	   collisions only, but let's treat all collision partners this way for now. */
	double beta = 8.629e-6;

	/* Also take some values from the bottom of page 4. Gamma = 2.15 at 10000 K and 1.58 at 1000
	   K. Do a linear interpolation. */
	double bigUpsilon10 = (T - 1000) / 9000 * 2.15 + (10000 - T) / 9000 * 1.58;

	Cvv(1, 0) = beta / sqrt(T) * bigUpsilon10 / the_gv(1);
	Cvv(0, 1) = Cvv(1, 0) * the_gv(1) / the_gv(0) *
	            exp(-(the_ev(1) - the_ev(0)) / Constant::BOLTZMAN / T);
	return Cvv;
}

Eigen::VectorXd TwoLevelHardcoded::alphav(double /* unused T */) const
{
	// There is no ion
	return Eigen::Vector2d::Zero();
}
