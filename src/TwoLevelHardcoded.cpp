#include "TwoLevelHardcoded.h"
#include "Constants.h"

TwoLevelHardcoded::TwoLevelHardcoded()
{
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

EVector TwoLevelHardcoded::ev() const { return the_ev; }

EVector TwoLevelHardcoded::gv() const { return the_gv; }

EMatrix TwoLevelHardcoded::avv() const
{
	EMatrix the_avv(2, 2);
	// clang-format off
	the_avv << 0, 0,
		   2.29e-6, 0;
	// clang-format on
	return the_avv;
}

EMatrix TwoLevelHardcoded::extraAvv() const { return EMatrix::Zero(2, 2); }

EMatrix TwoLevelHardcoded::cvv(double T, double /* unused ne */, double /* unused np */) const
{
	EMatrix Cvv = EMatrix::Zero(2, 2);

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

EVector TwoLevelHardcoded::sourcev(double T, double np, double ne) const
{
	// There is no ion
	return EVector::Zero(2);
}

EVector TwoLevelHardcoded::sinkv(double T, double n, double np, double ne) const
{
	return EVector::Zero(2);
}
