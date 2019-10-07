#include "SpecialFunctions.hpp"
#include "Constants.hpp"
#include "Error.hpp"
#include "TemplatedUtils.hpp"

#include <cmath>
#include <exception>

double SpecialFunctions::voigt(double a, double u)
{
	if (a <= 0.0)
		Error::runtime("Bad argument (a = " + std::to_string(a) +
		               ") . a should be positive.");

	double sigma = 1.0;
	double lg = 2.0 * M_SQRT2 * a;
	double xx = M_SQRT2 * u;
	int r = 5; // ask 1e-5 precision

	double x, y, k;
	x = xx / sigma / 1.41421356;
	y = lg / 2 / sigma / 1.41421356;
	double r0, r1;
	r0 = 1.51 * std::exp(1.144 * r);
	r1 = 1.60 * std::exp(0.554 * r);

	const double rrtpi = 0.56418958;
	double y0, y0py0, y0q;
	y0 = 1.5;
	y0py0 = y0 + y0;
	y0q = y0 * y0;
	double c[6] = {1.0117281,   -0.75197147,    0.012557727,
	               0.010022008, -0.00024206814, 0.00000050084806};
	double s[6] = {1.393237,     0.23115241,     -0.15535147,
	               0.0062183662, 0.000091908299, -0.00000062752596};
	double t[6] = {0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249};

	int j;
	int rg1, rg2, rg3;
	double abx, xq, yq, yrrtpi;
	double xlim0, xlim1, xlim2, xlim3, xlim4;
	double a0 = 0, d0 = 0, d2 = 0, e0 = 0, e2 = 0, e4 = 0, h0 = 0, h2 = 0, h4 = 0, h6 = 0;
	double p0 = 0, p2 = 0, p4 = 0, p6 = 0, p8 = 0, z0 = 0, z2 = 0, z4 = 0, z6 = 0, z8 = 0;
	double xp[6], xm[6], yp[6], ym[6];
	double mq[6], pq[6], mf[6], pf[6];
	double d, yf, ypy0, ypy0q;

	rg1 = 1;
	rg2 = 1;
	rg3 = 1;
	yq = y * y;
	yrrtpi = y * rrtpi;
	xlim0 = r0 - y;
	xlim1 = r1 - y;
	xlim3 = 3.097 * y - 0.45;
	xlim2 = 6.8 - y;
	xlim4 = 18.1 * y + 1.65;
	if (y <= 1e-6)
	{
		xlim1 = xlim0;
		xlim2 = xlim0;
	}
	abx = fabs(x);
	xq = abx * abx;
	if (abx > xlim0)
	{
		k = yrrtpi / (xq + yq);
	}
	else if (abx > xlim1)
	{
		if (rg1 != 0)
		{
			rg1 = 0;
			a0 = yq + 0.5;
			d0 = a0 * a0;
			d2 = yq + yq - 1.0;
		}
		d = rrtpi / (d0 + xq * (d2 + xq));
		k = d * y * (a0 + xq);
	}
	else if (abx > xlim2)
	{
		if (rg2 != 0)
		{
			rg2 = 0;
			h0 = 0.5625 + yq * (4.5 + yq * (10.5 + yq * (6.0 + yq)));
			h2 = -4.5 + yq * (9.0 + yq * (6.0 + yq * 4.0));
			h4 = 10.5 - yq * (6.0 - yq * 6.0);
			h6 = -6.0 + yq * 4.0;
			e0 = 1.875 + yq * (8.25 + yq * (5.5 + yq));
			e2 = 5.25 + yq * (1.0 + yq * 3.0);
			e4 = 0.75 * h6;
		}
		d = rrtpi / (h0 + xq * (h2 + xq * (h4 + xq * (h6 + xq))));
		k = d * y * (e0 + xq * (e2 + xq * (e4 + xq)));
	}
	else if (abx < xlim3)
	{
		if (rg3 != 0)
		{
			rg3 = 0;
			z0 = 272.1014 +
			     y * (1280.829 +
			          y * (2802.870 +
			               y * (3764.966 +
			                    y * (3447.629 +
			                         y * (2256.981 +
			                              y * (1074.409 +
			                                   y * (369.1989 +
			                                        y * (88.26741 +
			                                             y * (13.39880 + y)))))))));
			z2 = 211.678 +
			     y * (902.3066 +
			          y * (1758.336 +
			               y * (2037.310 +
			                    y * (1549.675 +
			                         y * (793.4273 +
			                              y * (266.2987 +
			                                   y * (53.59518 + y * 5.0)))))));
			z4 = 78.86585 +
			     y * (308.1852 +
			          y * (497.3014 +
			               y * (479.2576 +
			                    y * (269.2916 + y * (80.39278 + y * 10.0)))));
			z6 = 22.03523 +
			     y * (55.02933 + y * (92.75679 + y * (53.59518 + y * 10.0)));
			z8 = 1.496460 + y * (13.39880 + y * 5.0);
			p0 = 153.5168 +
			     y * (549.3954 +
			          y * (919.4955 +
			               y * (946.8970 +
			                    y * (662.8097 +
			                         y * (328.2151 +
			                              y * (115.3772 +
			                                   y * (27.93941 +
			                                        y * (4.264678 +
			                                             y * 0.3183291))))))));
			p2 = -34.16955 +
			     y * (-1.322256 +
			          y * (124.5975 +
			               y * (189.7730 +
			                    y * (139.4665 +
			                         y * (56.81652 +
			                              y * (12.79458 + y * 1.2733163))))));
			p4 = 2.584042 +
			     y * (10.46332 +
			          y * (24.01655 +
			               y * (29.81482 + y * (12.79568 + y * 1.9099744))));
			p6 = -0.07272979 + y * (0.9377051 + y * (4.266322 + y * 1.273316));
			p8 = 0.0005480304 + y * 0.3183291;
		}
		d = 1.7724538 / (z0 + xq * (z2 + xq * (z4 + xq * (z6 + xq * (z8 + xq)))));
		k = d * (p0 + xq * (p2 + xq * (p4 + xq * (p6 + xq * p8))));
	}
	else
	{
		ypy0 = y + y0;
		ypy0q = ypy0 * ypy0;
		k = 0.0;
		for (j = 0; j <= 5; j++)
		{
			d = x - t[j];
			mq[j] = d * d;
			mf[j] = 1.0 / (mq[j] + ypy0q);
			xm[j] = mf[j] * d;
			ym[j] = mf[j] * ypy0;
			d = x + t[j];
			pq[j] = d * d;
			pf[j] = 1.0 / (pq[j] + ypy0q);
			xp[j] = pf[j] * d;
			yp[j] = pf[j] * ypy0;
		}
		if (abx <= xlim4)
		{
			for (j = 0; j <= 5; j++)
			{
				k = k + c[j] * (ym[j] + yp[j]) - s[j] * (xm[j] - xp[j]);
			}
		}
		else
		{
			yf = y + y0py0;
			for (j = 0; j <= 5; j++)
			{
				k = k +
				    (c[j] * (mq[j] * mf[j] - y0 * ym[j]) + s[j] * yf * xm[j]) /
				                    (mq[j] + y0q) +
				    (c[j] * (pq[j] * pf[j] - y0 * yp[j]) - s[j] * yf * xp[j]) /
				                    (pq[j] + y0q);
			}
			k = y * k + std::exp(-xq);
		}
	}
	return k;
}

double SpecialFunctions::maxwellBoltzman(double v, double T, double m)
{
	double twokT = 2 * Constant::BOLTZMAN * T;
	double v2 = v * v;
	return Constant::FPI * std::pow(m / Constant::PI / twokT, 1.5) * v2 *
	       std::exp(-m * v2 / twokT);
}

double SpecialFunctions::planck(double nu, double T)
{
	constexpr double twoPlanckOverCsquare =
	                2 * Constant::PLANCK / Constant::LIGHT / Constant::LIGHT;
	return twoPlanckOverCsquare * nu * nu * nu /
	       expm1(Constant::PLANCK * nu / Constant::BOLTZMAN / T);
}

double SpecialFunctions::lorentz(double x, double gamma)
{
	return gamma / (x * x + gamma * gamma) / Constant::PI;
}

double SpecialFunctions::inverse_lorentz(double l, double gamma)
{
	return std::sqrt(gamma / Constant::PI / l - gamma * gamma);
}

double SpecialFunctions::lorentz_percentile(double p , double gamma)
{
	return gamma * std::tan(Constant::PI * (p - .5));
}

double SpecialFunctions::gauss(double x, double sigma)
{
	double x_sigma = x / sigma;
	return exp(-x_sigma * x_sigma / 2) / sigma / Constant::SQRT2PI;
}

double SpecialFunctions::inverse_gauss(double g, double sigma)
{
	return sigma * std::sqrt(-2 * std::log(Constant::SQRT2PI * sigma * g));
}

SpecialFunctions::LookupTable2D::LookupTable2D(std::function<double(double, double)> ff)
                : _ff{ff}
{
}

void SpecialFunctions::LookupTable2D::setup(const Array& xv, const Array& yv)

{
	_xv = xv;
	_yv = yv;
	size_t nx = xv.size();
	size_t ny = yv.size();

	_xMin = _xv[0];
	_xMax = _xv[nx - 1];
	_yMin = _yv[0];
	_yMax = _yv[ny - 1];

	// Tabulate the function on the given grids
	_fvv.resize(nx, ny);
	for (size_t x = 0; x < nx; x++)
		for (size_t y = 0; y < ny; y++)
			_fvv(x, y) = _ff(xv[x], yv[y]);
}

double SpecialFunctions::LookupTable2D::operator()(double x, double y) const
{
	// If out of the tabulated range, calculate explicitly. TODO: option to provide an
	// approximate "emergency" function that should be used when values are needed outside
	// of the tabulated range.
	if (x < _xMin || x > _xMax || y < _yMin || y > _yMax)
	{
		std::cout << "Lookup table out of range: << " << x << ", " << y << std::endl;
		return _ff(x, y);
	}

	size_t iRight = TemplatedUtils::index(x, _xv);
	if (iRight == 0)
		iRight = 1;
	size_t iLeft = iRight - 1;

	size_t iUpper = TemplatedUtils::index(y, _yv);
	if (iUpper == 0)
		iUpper = 1;
	size_t iLower = iUpper - 1;

	return TemplatedUtils::interpolateRectangular(
	                x, y, _xv[iLeft], _xv[iRight], _yv[iLower], _yv[iUpper],
	                _fvv(iLeft, iLower), _fvv(iRight, iLower), _fvv(iLeft, iUpper),
	                _fvv(iRight, iUpper));
}
