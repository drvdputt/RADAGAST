#include "TwoLevel.h"
#include "Constants.h"
#include "NumUtils.h"
#include "SpecialFunctions.h"

#include <exception>
#include <vector>

namespace {
	const double SQRT2PI = 2.50662827463;
	const double csquare_twoh = Constant::LIGHT * Constant::LIGHT / 2. / Constant::PLANCK;
}

TwoLevel::TwoLevel(const std::vector<double>& wavelength)
	: _wavelengthv(wavelength), _n(0), _nc(0), _T(0)
{
//	// toy model of CII 158 um from https://www.astro.umd.edu/~jph/N-level.pdf bottom of page 4
//	// this makes collisional deexcitation important for densities > 20-50 / cm3
//	// Level and transition information
//	_Ev << 0, .00786 / Constant::ERG_EV;
//	_gv << 1, 1;
//	_nv << _n, 0;

	// The only spontaneous transition is A_10
//	_Avv << 0, 0,
//		   2.29e-6, 0;

	// Lyman alpha (1s and 2p) model
	// Energy in cm^-1, so multiply with hc
	_Ev << 0., 82258.9191133 * Constant::PLANCKLIGHT;
	_gv << 1, 3;
	_nv << _n, 0;

	_Avv << 0, 0,
			6.2649e+08, 0;
}

void TwoLevel::doLevels(double n, double nc, double T, const vector<double>& isrf, const vector<double>& source, const vector<double>& sink)
{
	if (isrf.size() != _wavelengthv.size()) throw range_error("Given ISRF and wavelength vectors do not have the same size");
	if (source.size() != 2 || sink.size() != 2) throw range_error("Source and/or sink term vector(s) of wrong size");

	_n = n;
	_nc = nc;
	_T = T;

	// shortcut, in case of full ionization
	if (_n <= 0.)
	{
		_nv = Eigen::Vector2d::Zero();
	}
	else
	{
		// Calculate BijPij (needs to be redone at each temperature because the line profile can change)
		prepareAbsorptionMatrix(isrf);

		// Calculate Cij
		prepareCollisionMatrix();

		cout << "Aij" << endl << _Avv << endl << endl;
		cout << "BPij" << endl << _BPvv << endl << endl;
		cout << "Cij" << endl << _Cvv << endl << endl;

		// Calculate Fij and bi and solve F.n = b
		solveRateEquations(
				Eigen::Vector2d(source.data()),
				Eigen::Vector2d(sink.data()),
				0);
	}
}

double TwoLevel::bolometricEmission(size_t upper, size_t lower) const
{
	return (_Ev(upper) - _Ev(lower)) / Constant::FPI * _nv(upper) * _Avv(1,0);
}

std::vector<double> TwoLevel::calculateEmission() const
{
	// There is only one line for now
	double lineIntensity = bolometricEmission(1, 0);
	Eigen::ArrayXd result = lineIntensity * lineProfile(1, 0);
	std::vector<double> resultv(result.data(), result.data() + result.size());
	return resultv;
}

std::vector<double> TwoLevel::calculateOpacity() const
{
	double nu_ij = (_Ev(1) - _Ev(0)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT*Constant::LIGHT / 8. / Constant::PI / nu_ij/nu_ij * _Avv(1,0);
	double densityFactor = _nv(0)*_gv(1)/_gv(0) - _nv(1);
	Eigen::ArrayXd result = constantFactor * densityFactor * lineProfile(1, 0);
	return std::vector<double>(result.data(), result.data() + result.size());
}

Eigen::ArrayXd TwoLevel::lineProfile(size_t upper, size_t lower) const
{
	double nu0 = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;

	double decayRate = _Avv(upper, lower) + _Cvv(upper, lower) // decay rate of top level
			+ _Cvv(lower, upper); // decay rate of bottom level

	double thermalVelocity = std::sqrt(Constant::BOLTZMAN * _T / Constant::HMASS_CGS);

	// Half the FWHM of the Lorentz
	double halfWidth = decayRate / Constant::FPI;

	// The standard deviation in frequency units. It is about half of the FWHM for a Gaussian
	double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;
	double one_sqrt2sigma = M_SQRT1_2 / sigma_nu;

	Eigen::ArrayXd profile(_wavelengthv.size());
	for (size_t n = 0; n < _wavelengthv.size(); n++)
	{
		double deltaNu = Constant::LIGHT / _wavelengthv[n] - nu0;

		// When we are very far in the tail, just skip the calculation
		if (false)/*(std::abs(deltaNu) > 20 * sumWidths)*/
		{
			profile(n) = 0.;
		}
		else
		{
			profile(n) = SpecialFunctions::voigt(one_sqrt2sigma * halfWidth, one_sqrt2sigma * deltaNu) / SQRT2PI / sigma_nu;
			// phi_lambda dlambda = -phi_nu dnu
			// phi_lambda = -phi_nu * d(c / lambda)/dlambda = phi_nu * c / lambda^2
			profile(n) *= Constant::LIGHT / _wavelengthv[n] / _wavelengthv[n];
		}
	}
	cout << "line profile norm = " << NumUtils::integrate<double>(_wavelengthv, std::vector<double>(profile.data(), profile.data() + profile.size())) << endl;
	return profile;
}

void TwoLevel::prepareAbsorptionMatrix(const std::vector<double>& isrf)
{
	// Calculate product of Bij and Pij = integral(phi * I_lambda)
	for (int i = 0; i < _BPvv.rows(); i++)
	{
		// Go over the lower triangle, without diagonal
		for(int j = 0; j < i; j++)
		{
			// Calculate Pij for the lower triangle (= stimulated emission)
			Eigen::ArrayXd phi = lineProfile(i, j);
			std::vector<double> integrand;
			integrand.reserve(isrf.size());
			for (size_t n = 0; n < isrf.size(); n++)
			{
				// since our formula for B assumes a specific intensity in frequency units as measure of the radiation field,
				// we need to compute integral I_nu(lambda) phi_lambda(lambda) dlambda
				// I_nu(lambda) = lambda * lambda * I_lambda(lambda) / c
				// = lambda * lambda * u_lambda(lambda) / 4pi
				// = lambda * lambda * isrf / 4pi
				integrand.push_back(_wavelengthv[n]*_wavelengthv[n] * phi(n) * isrf[n]);
			}
			_BPvv(i, j) = NumUtils::integrate<double>(_wavelengthv, integrand) / Constant::FPI;

			// Multiply by Bij in terms of Aij, valid for i > j
			double nu_ij = (_Ev(i) - _Ev(j)) / Constant::PLANCK;
			_BPvv(i, j) *= csquare_twoh / nu_ij / nu_ij / nu_ij * _Avv(i,j);

			// Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij = Pji
			_BPvv(j, i) = _gv(i) / _gv(j) * _BPvv(i,j);
		}
		// Set the diagonal to zero
		_BPvv(i, i) = 0;
	}
}

void TwoLevel::prepareCollisionMatrix()
{
//	// Need separate contributions for number of protons and electrons
//	// Toy implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
//	// is actually for electron collisions only, but let's treat all collision partners this way for now
//	double beta = 8.629e-6;
//
//	// also take some values from the bottom of page 4
//	// Gamma = 2.15 at 10000 K and 1.58 at 1000 K
//	double bigGamma10 = (_T - 1000) / 9000 * 2.15 + (10000 - _T) / 9000 * 1.58;
//
//	double C10 = beta / sqrt(_T) * bigGamma10 / _gv(1) * _nc;
//	double C01 = C10 * _gv(1) / _gv(0) * exp(-(_Ev(1) - _Ev(0)) / Constant::BOLTZMAN / _T);
//	_Cvv << 0, C01,
//			C10, 0;
	_Cvv << 0,0,0,0;
}

void TwoLevel::solveRateEquations(Eigen::Vector2d sourceTerm, Eigen::Vector2d sinkTerm, size_t chooseConsvEq)
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	Eigen::MatrixXd Mvv(_Avv.transpose() + _BPvv.transpose() + _Cvv.transpose());

	// See equation for Fij (37) in document
	// = subtract departure rate from level i to all other levels
	Eigen::MatrixXd departureDiagonal = Mvv.colwise().sum().asDiagonal();
	cout << "Mij" << endl << Mvv << endl << endl;
	cout << "departure" << endl << departureDiagonal << endl << endl;
	Mvv -= departureDiagonal;
	Mvv -= sinkTerm.asDiagonal();
	Eigen::VectorXd b(-sourceTerm);

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = Eigen::VectorXd::Ones(Mvv.cols());
	b(chooseConsvEq) = _n;

	std::cout << "System to solve:\n" << Mvv << " * nv\n=\n" << b << endl << endl;

	_nv = Mvv.colPivHouseholderQr().solve(b);
	cout << "nv" << endl << _nv << endl;
	double relative_error = (Mvv*_nv - b).norm() / b.norm(); // norm() is L2 norm
	cout << "The relative error is: " << relative_error << endl;

	// use explicit formula when this row has very large coefficients (only works if chooseConsvEq = 0)
	// otherwise there can be problems when subtracting the doubles from each other
	// not sure how I will solve this for more general systems
	if (abs(max( Mvv(1, 0), Mvv(1, 1) )) > 1e15 )
	{
		cout << "Overriding with explicit solution" << endl;
		_nv(0) = (-Mvv(1, 1)*_n - sinkTerm(1)) / (-Mvv(1, 1) + Mvv(1, 0));
		_nv(1) = _n - _nv(0);
		cout << "nv" << endl << _nv << endl;
		relative_error = (Mvv*_nv - b).norm() / b.norm(); // norm() is L2 norm
		cout << "The relative error is: " << relative_error << endl;
	}
	cout << "from matrix equation nu / nl: " << _nv(1) / _nv(0) << endl;
	cout << "Directly from coefficients (no sources) nu / nl: " << (_Cvv(0,1) + _BPvv(0,1)) / (_Avv(1,0) + _Cvv(1,0) + _BPvv(1,0)) << endl;
}




