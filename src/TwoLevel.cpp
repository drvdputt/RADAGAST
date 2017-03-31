#include "TwoLevel.h"
#include "Constants.h"
#include "Sanity.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include <algorithm>
#include <exception>
#include <vector>
#include "global.h"

using namespace std;

#define SOLUTION_TOLERANCE 0.1

TwoLevel::TwoLevel(const Array& frequencyv) : _frequencyv(frequencyv)
{
	//	// toy model of CII 158 um from https://www.astro.umd.edu/~jph/N-level.pdf bottom of
	// page 4
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

	_Avv << 0, 0, 6.2649e+08, 0;
}

void TwoLevel::solveBalance(double n, double ne, double np, double T,
                            const Array& specificIntensity, const Array& source, const Array& sink)
{
	if (specificIntensity.size() != _frequencyv.size())
		throw range_error("Given ISRF and wavelength vectors do not have the same size");
	if (source.size() != 2 || sink.size() != 2)
		throw range_error("Source and/or sink term vector(s) of wrong size");

	_n = n;
	_ne = ne;
	_np = np;
	_T = T;

	if (_n > 0)
	{
		// Calculate Cij (needs to happen bfore the first call to the line profile. Maybe
		// the dependencies need to be made clearer by adopting a more functional
		// programming-like style)
		prepareCollisionMatrix();

#ifdef REPORT_LINE_QUALITY
		double norm = TemplatedUtils::integrate<double>(_frequencyv, lineProfile(1, 0));
		DEBUG("line profile norm = " << norm << endl);
#endif

		// Calculate BijPij (needs to be redone at each temperature because the line profile
		// can change) Also needs the Cij to calculate collisional broadening
		prepareAbsorptionMatrix(specificIntensity);

#ifdef PRINT_MATRICES
		DEBUG("Aij" << endl << _Avv << endl << endl);
		DEBUG("BPij" << endl << _BPvv << endl << endl);
		DEBUG("Cij" << endl << _Cvv << endl << endl);
#endif
		// Calculate Fij and bi and solve F.n = b
		solveRateEquations(Eigen::Map<const Eigen::Vector2d>(&source[0]),
		                   Eigen::Map<const Eigen::Vector2d>(&sink[0]), 0);
	}
	else
	{
		_nv = Eigen::Vector2d::Zero();
	}
}

Array TwoLevel::emissivityv() const { return lineIntensityFactor(1, 0) * lineProfile(1, 0); }

Array TwoLevel::opacityv() const { return lineOpacityFactor(1, 0) * lineProfile(1, 0); }

Array TwoLevel::scatteringOpacityv() const { return opacityv() * lineDecayFraction(1, 0); }

double TwoLevel::lineIntensityFactor(size_t upper, size_t lower) const
{
	return (_Ev(upper) - _Ev(lower)) / Constant::FPI * _nv(upper) * _Avv(1, 0);
}

double TwoLevel::lineOpacityFactor(size_t upper, size_t lower) const
{
	double nu_ij = (_Ev(1) - _Ev(0)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij /
	                        nu_ij * _Avv(1, 0);
	double densityFactor = _nv(0) * _gv(1) / _gv(0) - _nv(1);
	return constantFactor * densityFactor;
}

Array TwoLevel::lineProfile(size_t upper, size_t lower) const
{
	double nu0 = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;

	double decayRate = _Avv(upper, lower) + _Cvv(upper, lower) // decay rate of top level
	                   + _Cvv(lower, upper);                   // decay rate of bottom level
	// (stimulated emission doesn't count)

	double thermalVelocity = sqrt(Constant::BOLTZMAN * _T / Constant::HMASS_CGS);

	// Half the FWHM of the Lorentz
	double halfWidth = decayRate / Constant::FPI;

	// The standard deviation in frequency units. It is about half of the FWHM for a Gaussian
	double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;
	double one_sqrt2sigma = M_SQRT1_2 / sigma_nu;

	Array profile(_frequencyv.size());
	for (size_t n = 0; n < _frequencyv.size(); n++)
	{
		double deltaNu = _frequencyv[n] - nu0;
		profile[n] = SpecialFunctions::voigt(one_sqrt2sigma * halfWidth,
		                                     one_sqrt2sigma * deltaNu) /
		             Constant::SQRT2PI / sigma_nu;
	}
	return profile;
}

double TwoLevel::lineDecayFraction(size_t upper, size_t lower) const
{
	return _Avv(upper, lower) /
	       (_Avv.row(upper).sum() + _BPvv.row(upper).sum() + _Cvv.row(upper).sum());
}

void TwoLevel::prepareAbsorptionMatrix(const Array& specificIntensityv)
{
	// Calculate product of Bij and Pij = integral(phi * I_nu)
	for (int i = 0; i < _BPvv.rows(); i++)
	{
		// Go over the lower triangle, without diagonal
		for (int j = 0; j < i; j++)
		{
			// Calculate Pij for the lower triangle (= stimulated emission)
			_BPvv(i, j) = TemplatedUtils::integrate<double, Array, Array>(
			                _frequencyv, lineProfile(i, j) * specificIntensityv);

			// Multiply by Bij in terms of Aij, valid for i > j
			double nu_ij = (_Ev(i) - _Ev(j)) / Constant::PLANCK;
			_BPvv(i, j) *= Constant::CSQUARE_TWOPLANCK / nu_ij / nu_ij / nu_ij *
			               _Avv(i, j);

			// Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij =
			// Pji
			_BPvv(j, i) = _gv(i) / _gv(j) * _BPvv(i, j);
		}
		// Set the diagonal to zero
		_BPvv(i, i) = 0;
	}
}

void TwoLevel::prepareCollisionMatrix()
{
	//	// Need separate contributions for number of protons and electrons
	//	// Toy implementation below, inspired by https://www.astro.umd.edu/~jph/N-level.pdf
	//	// is actually for electron collisions only, but let's treat all collision partners
	// this way for now 	double beta = 8.629e-6;
	//
	//	// also take some values from the bottom of page 4
	//	// Gamma = 2.15 at 10000 K and 1.58 at 1000 K
	//	double bigGamma10 = (_T - 1000) / 9000 * 2.15 + (10000 - _T) / 9000 * 1.58;
	//
	//	// data from 2002-anderson
	// CORRECTION: need i = 3 is for the 2p state
	vector<double> electronTemperaturesEV = {.5, 1., 3., 5., 10., 15., 20., 25.};

	//	// i = 2
	//	vector<double> effectiveCollisionStrv = {2.6e-1,  2.96e-1, 3.26e-1, 3.39e-1,
	//	                                         3.73e-1, 4.06e-1, 4.36e-1, 4.61e-1};

	// i = 3
	vector<double> effectiveCollisionStrv = {4.29e-1, 5.29e-01, 8.53e-01, 1.15e00,
	                                         1.81e00, 2.35e00,  2.81e00,  3.20e00};

	double currentT_eV = Constant::BOLTZMAN * _T * Constant::ERG_EV;

	// naively inter- and extrapolate this data linearly, just to have some form of collision
	// coefficient (provides numerical stability)
	size_t iRight = TemplatedUtils::index(currentT_eV, electronTemperaturesEV);
	if (iRight == 0)
		iRight = 1;
	else if (iRight == effectiveCollisionStrv.size())
		iRight -= 1;

	// Can be larger than one or lower than 0 when the given temperature is outside of the data
	// range. In that case, a linear extrapolation is the result.
	double wRight = (currentT_eV - electronTemperaturesEV[iRight - 1]) /
	                (electronTemperaturesEV[iRight] - electronTemperaturesEV[iRight - 1]);

	double bigGamma10 = (1. - wRight) * effectiveCollisionStrv[iRight - 1] +
	                    wRight * effectiveCollisionStrv[iRight];

	const double I_H_eV = 13.6058;
	// n-changing collisions are dominated by electrons, so multiply with _ne
	double C10 = 2.1716e-8 / _gv(1) * sqrt(I_H_eV / currentT_eV) * bigGamma10 * _ne;

	double C01 = C10 * _gv(1) / _gv(0) * exp(-(_Ev(1) - _Ev(0)) / Constant::BOLTZMAN / _T);
	_Cvv << 0, C01, C10, 0;
}

void TwoLevel::solveRateEquations(Eigen::Vector2d sourceTerm, Eigen::Vector2d sinkTerm,
                                  size_t chooseConsvEq)
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	Eigen::MatrixXd Mvv(_Avv.transpose() + _BPvv.transpose() + _Cvv.transpose());

	// See equation for Fij (37) in document
	// subtract departure rate from level i to all other levels
	Eigen::MatrixXd departureDiagonal = Mvv.colwise().sum().asDiagonal();
	Mvv -= departureDiagonal;
	Mvv -= sinkTerm.asDiagonal();
	Eigen::VectorXd b(-sourceTerm);

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = Eigen::VectorXd::Ones(Mvv.cols());
	b(chooseConsvEq) = _n;

#ifdef PRINT_MATRICES
	DEBUG("System to solve:\n" << Mvv << " * nv\n=\n" << b << endl << endl);
#endif

	// Call the linear solver
	_nv = Mvv.colPivHouseholderQr().solve(b);

	DEBUG("nv" << endl << _nv << endl);
	double err0, err1;
	Eigen::Vector2d diff = Mvv * _nv - b;
	err0 = diff(0) / b(0);
	err1 = diff(1) / b(1);
	DEBUG("The relative errors are: " << err0 << " and " << err1 << endl);

	// use explicit formula when the solution is clearly wrong
	// not sure how I will solve this for more general systems
	if (abs(err0) > SOLUTION_TOLERANCE || abs(err1) > SOLUTION_TOLERANCE)
	{
		//		double nv0 = (-Mvv(1, 1) * _n - b(1) / (-Mvv(1, 1) + Mvv(1, 0));
		//		double nv1 = _n - nv0;
		double nv1 = (Mvv(1, 0) * _n - b(1)) / (Mvv(1, 0) - Mvv(1, 1));
		double nv0 = _n - nv1;
#ifdef REPORT_OVERRIDE
		Sanity::reportOverridden("nv(0)", _nv(0), nv0, "numerical errors where too high.");
		Sanity::reportOverridden("nv(1)", _nv(1), nv1, "numerical errors where too high.");
#endif
		_nv(0) = nv0;
		_nv(1) = nv1;
		diff = Mvv * _nv - b;
		err0 = diff(0) / b(0);
		err1 = diff(1) / b(1);
		DEBUG("The relative errors are: " << err0 << " and " << err1 << endl);
	}
	//	DEBUG("from matrix equation nu / nl: " << _nv(1) / _nv(0) << endl);
	//	DEBUG("Directly from coefficients (no sources) nu / nl: "
	//			<< (_Cvv(0, 1) + _BPvv(0, 1)) / (_Avv(1, 0) + _Cvv(1, 0) + _BPvv(1,
	// 0))
	//<< endl);
}
