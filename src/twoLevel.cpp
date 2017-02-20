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
	: _wavelength(wavelength), _n(0), _nc(0), _T(0)
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

void TwoLevel::doLevels(double n, double nc, double T, const std::vector<double>& isrf, double recombinationRate)
{
	if (isrf.size() != _wavelength.size()) throw std::range_error("Given ISRF and wavelength vectors do not have the same size");

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
		// The ionization rate calculation makes no distinction between the levels.
		// Therefore, the sink term is the same for each level. Moreover, total source = total sink
		// so we want sink*n0 + sink*n1 = source => sink = source / n because n0/n + n1/n = 1
		solveRateEquations(Eigen::Vector2d(0.1*recombinationRate, 0.9*recombinationRate),
				Eigen::Vector2d(recombinationRate/_n, recombinationRate/_n), 0);
	}
}

double TwoLevel::bolometricEmission(size_t upper, size_t lower) const
{
	return (_Ev(upper) - _Ev(lower)) * Constant::FPI * _nv(upper) * _Avv(1,0);
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
	double constantFactor = Constant::LIGHT*Constant::LIGHT / 8. / Constant::PI / nu_ij*nu_ij * _Avv(1,0);
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
			profile(n) = SpecialFunctions::voigt(one_sqrt2sigma * halfWidth, one_sqrt2sigma * deltaNu) / SQRT2PI / sigma_nu;
			// phi_lambda dlambda = -phi_nu dnu
			// phi_lambda = -phi_nu * d(c / lambda)/dlambda = phi_nu * c / lambda^2
			profile(n) *= Constant::LIGHT / _wavelength[n] / _wavelength[n];
		}
	}
	cout << "line profile norm = " << NumUtils::integrate<double>(_wavelength, std::vector<double>(profile.data(), profile.data() + profile.size())) << endl;
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
			for (size_t n = 0; n < isrf.size(); n++) integrand.push_back(phi(n) * isrf[n]);
			// Make sure that isrf is in specific intensity units
			_BPvv(i, j) = Constant::LIGHT / Constant::FPI * NumUtils::integrate<double>(_wavelength, integrand);

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
	// Need separate contributions for number of protons and electrons
	// Toy implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
	// is actually for electron collisions only, but let's treat all collision partners this way for now
	double beta = 8.629e-6;

	// also take some values from the bottom of page 4
	// Gamma = 2.15 at 10000 K and 1.58 at 1000 K
	double bigGamma10 = (_T - 1000) / 9000 * 2.15 + (10000 - _T) / 9000 * 1.58;

	double C10 = beta / sqrt(_T) * bigGamma10 / _gv(1) * _nc;
	double C01 = C10 * _gv(1) / _gv(0) * exp(-(_Ev(1) - _Ev(0)) / Constant::BOLTZMAN / _T);
	_Cvv << 0, C01,
			C10, 0;
}

void TwoLevel::solveRateEquations(Eigen::Vector2d sourceTerm, Eigen::Vector2d sinkTerm, size_t chooseConsvEq)
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	Eigen::MatrixXd Mvv(_Avv.transpose() + _BPvv.transpose() + _Cvv.transpose());

	// See equation for Fij (37) in document
	// = subtract departure rate from level i to all other levels
	Eigen::MatrixXd departureDiagonal = Mvv.colwise().sum().asDiagonal();
	cout << "departure" << endl << departureDiagonal << endl << endl;
	Mvv -= Mvv.colwise().sum().asDiagonal();
	Mvv -= sinkTerm.asDiagonal();
	Eigen::VectorXd b(-sourceTerm);

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = Eigen::VectorXd::Ones(Mvv.cols());
	b(chooseConsvEq) = _n;

	std::cout << "System to solve:\n" << Mvv << "\n" << b << std::endl;

	_nv = Mvv.colPivHouseholderQr().solve(b);

	cout << "nv" << endl << _nv << endl;
	cout << "from matrix equation nu / nl: " << _nv(1) / _nv(0) << endl;

	//cout << "from explicit formula nu / nl: " << _gv(1) / _gv(0) * exp((_Ev(0) - _Ev(1)) / _T / Constant::BOLTZMAN) / (1 + _Avv(1,0) / _Cvv(1,0)) << endl;
	cout << "from explicit formula nu / nl: " << (_Cvv(0,1) + _BPvv(0,1)) / (_Avv(1,0) + _Cvv(1,0) + _BPvv(1,0)) << endl;
}





