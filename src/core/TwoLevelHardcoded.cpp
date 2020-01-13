#include "TwoLevelHardcoded.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"

TwoLevelHardcoded::TwoLevelHardcoded() : LevelCoefficients(12 * Constant::HMASS_CGS)
{
    // Level and transition information
    //	_Ev << 0, .00786 / Constant::ERG_EV;
    //	_gv << 1, 1;
    //	_nv << _n, 0;

    // The only spontaneous transition is A_10
    //	_Avv << 0, 0,
    //		   2.29e-6, 0;
    EVector the_ev(2);
    EVector the_gv(2);
    the_ev << 0, .00786 / Constant::ERG_EV;
    the_gv << 1, 1;

    EMatrix the_avv(2, 2);
    // clang-format off
	the_avv << 0, 0,
		   2.29e-6, 0;
    // clang-format on
    setConstants(the_ev, the_gv, the_avv, EMatrix::Zero(2, 2));
}

EMatrix TwoLevelHardcoded::cvv(const CollisionParameters& cp) const
{
    EMatrix Cvv = EMatrix::Zero(2, 2);
    double ne = cp._sv.ne();
    if (ne > 0)
    {
        /* Need separate contributions for number of protons and electrons. Toy
		   implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
		   is actually for electron collisions only. */
        double beta = 8.629e-6;

        /* Also take some values from the bottom of page 4. Gamma = 2.15 at 10000 K and
		   1.58 at 1000 K. Do a linear interpolation. */
        double T = cp._t;
        double bigUpsilon10 = (T - 1000) / 9000 * 2.15 + (10000 - T) / 9000 * 1.58;

        Cvv(1, 0) += beta / sqrt(T) * bigUpsilon10 / gv()(1) * ne;
        Cvv(0, 1) += Cvv(1, 0) * gv()(1) / gv()(0) * exp(-(ev()(1) - ev()(0)) / Constant::BOLTZMAN / T);
    }
    return Cvv;
}
