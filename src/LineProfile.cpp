#include "LineProfile.h"
#include "Constants.h"
#include "SpecialFunctions.h"

#include <cmath>

LineProfile::LineProfile(double center, double sigma_gauss, double halfWidth_lorenz)
                : _center{center}, _sigma_gauss{sigma_gauss}

{
	_one_sqrt2sigma = M_SQRT1_2 / _sigma_gauss;
	_a = halfWidth_lorenz * _one_sqrt2sigma;
}

double LineProfile::operator()(double nu) const
{
	double x = (nu - _center) * _one_sqrt2sigma;
	// Note that the normalization factor is 1 / sqrt(2 pi sigma)
	return SpecialFunctions::voigt(_a, x) / Constant::SQRT2PI / _sigma_gauss;
}
