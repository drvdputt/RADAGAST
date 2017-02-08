#include "TwoLevel.h"
#include "Constants.h"
#include "SpecialFunctions.h"

#include <vector>

namespace {
	const double SQRT2PI = 2.50662827463;
}

TwoLevel::TwoLevel(double n, const std::vector<double>& wavelength, const std::vector<double>& isrf)
	: _n(n), _wavelength(wavelength), _isrf(isrf), _T(0), _n0(n), _n1(0)
{
	// Level and transition information
	// The only spontaneous transition is A_10
	Aij << 0, 0,
		   1, 0;

	_E0 = 0;
	_E1 = 10. / Constant::ERG_EV;
	_g0 = 1;
	_g1 = 3;
}

void TwoLevel::doLevels(double T)
{
	_T = T;
	Eigen::Matrix2d Mij = constructRateMatrix();
	Eigen::Vector2d nv = solveRateEquations(Mij, Eigen::Vector2d::Zero(2));
	_n0 = nv(0);
	_n1 = nv(1);
}

double TwoLevel::bolometricEmission() const
{
	double deltaE = _E1 - _E0;
	return deltaE / Constant::FPI * _n1 * Aij(1,0);
}

std::vector<double> TwoLevel::calculateEmission() const
{
	// There is only one line for now
	double deltaE = _E1 - _E0;
	double lineIntensity = deltaE / Constant::FPI * _n1 * Aij(1,0);
	Eigen::ArrayXd result = lineIntensity * lineProfile(0);
	return std::vector<double>(result.data(), result.data() + result.size());
}

std::vector<double> TwoLevel::calculateOpacity() const
{
	double nu_ij = (_E1 - _E0) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT*Constant::LIGHT / 8. / Constant::PI / nu_ij*nu_ij * Aij(1,0);
	double densityFactor = _n0*_g1/_g0 - _n1;
	Eigen::ArrayXd result = constantFactor * densityFactor * lineProfile(0);
	return std::vector<double>(result.data(), result.data() + result.size());
}

Eigen::ArrayXd TwoLevel::lineProfile(size_t lineIndex) const
{
	double nu0 = (_E1 - _E0) / Constant::PLANCK;

	double decayRate = Aij(1,0) + Cij(1,0) // decay rate of top level
			+ Cij(0,1); // decay rate of bottom level

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

Eigen::Matrix2d TwoLevel::constructRateMatrix()
{
	return Eigen::Matrix2d();
}

Eigen::Vector2d TwoLevel::solveRateEquations(Eigen::Matrix2d Mij, Eigen::Vector2d alpha)
{
	return Eigen::Vector2d();
}





