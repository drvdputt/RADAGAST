#include "TwoLevel.h"
#include "Constants.h"
#include "NumUtils.h"
#include "SpecialFunctions.h"

#include <vector>

namespace {
	const double SQRT2PI = 2.50662827463;
}

TwoLevel::TwoLevel(double n, const std::vector<double>& wavelength, const std::vector<double>& isrf)
	: _n(n), _isrf(isrf), _wavelength(wavelength), _T(0)
{
	// toy model of CII 158 um from https://www.astro.umd.edu/~jph/N-level.pdf bottom of page 4
	// this makes collisional deexcitation important for densities > 20-50 / cm3
	// Level and transition information
	_Ev << 0, .00786 / Constant::ERG_EV;
	_gv << 1, 1;
	_nv << _n, 0;

	// The only spontaneous transition is A_10
	_Avv << 0, 0,
		   2.29e-6, 0;
}

void TwoLevel::doLevels(double T)
{
	_T = T;

	// Calculate Pij = integral phi * I_lambda
	for (size_t i = 0; i < _Pvv.rows(); i++)
	{
		for(size_t j = i + 1; j < _Pvv.cols(); j++)
		{
			Eigen::ArrayXd phi = lineProfile(i, j);
			std::vector<double> integrand;
			integrand.reserve(_isrf.size());
			for (size_t n = 0; n < integrand.size(); n++) integrand.push_back(phi(n) * _isrf[n]);
			_Pvv(i, j) = NumUtils::integrate<double>(_wavelength, integrand);

			// Copy to the lower triangle
			_Pvv(j, i) = _Pvv(i, j);
		}
		// Set the diagonal to zero
		_Pvv(i, i) = 0;
	}

	// Calculate Cij
	prepareCollisionMatrix();

	cout << "Aij" << endl << _Avv << endl;
	cout << "Pij" << endl << _Pvv << endl;
	cout << "Cij" << endl << _Cvv << endl;

	// Calculate Fij and bi and solve F.n = b
	solveRateEquations(Eigen::Vector2d::Zero(2), 0);
}

double TwoLevel::bolometricEmission() const
{
	return (_Ev(1) - _Ev(0)) * Constant::FPI * _nv(1) * _Avv(1,0);
}

std::vector<double> TwoLevel::calculateEmission() const
{
	// There is only one line for now
	double lineIntensity = (_Ev(1) - _Ev(0)) / Constant::FPI * _nv(1) * _Avv(1,0);
	Eigen::ArrayXd result = lineIntensity * lineProfile(0, 1);
	return std::vector<double>(result.data(), result.data() + result.size());
}

std::vector<double> TwoLevel::calculateOpacity() const
{
	double nu_ij = (_Ev(1) - _Ev(0)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT*Constant::LIGHT / 8. / Constant::PI / nu_ij*nu_ij * _Avv(1,0);
	double densityFactor = _nv(0)*_gv(1)/_gv(0) - _nv(0);
	Eigen::ArrayXd result = constantFactor * densityFactor * lineProfile(0, 1);
	return std::vector<double>(result.data(), result.data() + result.size());
}

Eigen::ArrayXd TwoLevel::lineProfile(size_t i, size_t j) const
{
	double nu0 = (_Ev(i) - _Ev(j)) / Constant::PLANCK;

	double decayRate = _Avv(1,0) + _Cvv(1,0) // decay rate of top level
			+ _Cvv(0,1); // decay rate of bottom level

	double thermalVelocity = std::sqrt(Constant::BOLTZMAN * _T / Constant::HMASS_CGS);

	// Half the FWHM of the Lorentz
	double halfWidth = decayRate / Constant::FPI;

	// The standard deviation in frequency units. It is about half of the FWHM for a Gaussian
	double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;

	// Get an order of magnitude guess of the total width: Gamma / 2 + sigma_nu
	double sumWidths = halfWidth + sigma_nu;

	Eigen::ArrayXd profile(_wavelength.size());
	for (size_t n = 0; n < _wavelength.size(); n++)
	{
		double deltaNu = Constant::LIGHT / _wavelength[n] - nu0;

		// When we are very far in the tail, just skip the calculation
		if (std::abs(deltaNu) > 20 * sumWidths) profile(n) = 0.;

//		// If natural+collisional effect dominates, ignore the Gaussian component.
//		if (halfWidth > 1e4 * sigma_nu)
//		{
//			return halfWidth / Constant::PI / (deltaNu*deltaNu + halfWidth*halfWidth); // Lorentz
//		}
//		// If thermal+turbulent effect dominates, ignore the Lorentzian component.
//		else if (decayRate < 1e-4 * sigma_nu)
//		{
//			return std::exp(-deltaNu*deltaNu / 2. / sigma_nu/sigma_nu) / SQRT2PI / sigma_nu; // Gauss
//		}
//		// When both components are important, use this approximation of the Voigt profile
		else
		{
			double one_sqrt2sigma = M_SQRT1_2 / sigma_nu;
			profile(n) = SpecialFunctions::voigt(one_sqrt2sigma * halfWidth, one_sqrt2sigma * deltaNu) / SQRT2PI / sigma_nu;
		}
	}
	return profile;
}

void TwoLevel::prepareCollisionMatrix()
{
	// Need separate contributions for number of protons and electrons
	// Toy implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
	// is actually for electron collisions only, but let's treat all collision partners this way for now
	double beta = 8.629e-6;

	// also take some values from the bottom of page 4
	// Gamma = 2.15 at 10000 K and 1.58 at 1000 K
	double bigGamma10 = (_T - 1000) / 9000 * 2.15 + (10000 - _T) / 9000 * 1.58;

	double C10 = beta / sqrt(_T) * bigGamma10 / _gv(1) * _n;
	double C01 = C10 * _gv(1) / _gv(0) * exp(-(_Ev(1) - _Ev(0)) / Constant::BOLTZMAN / _T);
	_Cvv << 0, C01,
			C10, 0;
}

void TwoLevel::solveRateEquations(Eigen::Vector2d ne_np_alpha, size_t chooseConsvEq)
{
	const double csquare_twoh = Constant::LIGHT * Constant::LIGHT / 2. / Constant::PLANCK;

	// Initialize Mij as Aji + Cji
	Eigen::MatrixXd Mvv(_Avv.transpose() + _Cvv.transpose());

	// Add BjiPji (loop over the upper triangle j > i of M(i,j) (= lower triangle of B(j,i)) (no diagonal elements, )
	for (size_t i = 0; i < Mvv.rows(); i++)
	{
		for(size_t j = i + 1; j < Mvv.cols(); j++)
		{
			double nu_ji = (_Ev(j) - _Ev(i)) / Constant::PLANCK;

			// Expression for Bji (firstindex > secondindex)
			double Bji = csquare_twoh / nu_ji / nu_ji / nu_ji * _Avv(j,i);
			Mvv(i,j) += Bji * _Pvv(j,i);

			// For Bij (firstindex < secondindex)
			Mvv(j,i) += _gv(j) / _gv(i) * Bji * _Pvv(i,j);
		}
	}

	// See equation for Fij (37) in document
	Mvv -= Mvv.colwise().sum().asDiagonal();
	Eigen::VectorXd b(-ne_np_alpha);

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = Eigen::VectorXd::Ones(Mvv.cols());
	b(chooseConsvEq) = _n;

	std::cout << "System to solve:\n" << Mvv << "\n" << b << std::endl;

	_nv = Mvv.colPivHouseholderQr().solve(b);

	cout << "nv" << endl << _nv << endl;
}





