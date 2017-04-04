#include "NLevel.h"
#include "Constants.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"
#include "global.h"
#include <iostream>

using namespace std;

namespace
{
// commonly used indices
const int index2p = 1;
const int index2s = 2;
}

NLevel::NLevel() : _Ev(_nLv), _gv(_nLv), _Avv(_nLv, _nLv), _extraAvv(_nLv, _nLv)
{
	// Using NIST data
	// (http://physics.nist.gov/cgi-bin/ASD/lines1.pl?unit=1&line_out=0&bibrefs=1&show_obs_wl=1&show_calc_wl=1&A_out=0&intens_out=1&allowed_out=1&forbid_out=1&conf_out=1&term_out=1&enrg_out=1&J_out=1&g_out=0&spectra=H%20I)
	// 1s, 2p, 2s
	// 3, 4, 5
	_Ev << 0., 82258.9191133, 82258.9543992821, 97492.304, 102823.904, 105291.657;
	// Energy in cm^-1, so multiply with hc
	_Ev *= Constant::PLANCKLIGHT;

	_gv << 1, 1, 3, 9, 16, 25;

	double A2p1 = 6.2649e8;

	double A2s1 = 2.495e-6;

	double A31 = 5.5751e7;
	double A32 = 4.4101e+07;
	double A32p = 3 * A32 / 4.;
	double A32s = A32 / 4.;

	double A41 = 1.2785e+07;
	double A42 = 8.4193e+06;
	double A42p = 3 * A42 / 4.;
	double A42s = A42 / 4.;
	double A43 = 8.9860e+06;

	double A51 = 4.1250e+06;
	double A52 = 2.5304e+06;
	double A52p = 3 * A52 / 4.;
	double A52s = A52 / 4.;
	double A53 = 2.2008e+06;
	double A54 = 2.6993e+06;

	// clang-format off
	_Avv << 0, 0, 0, 0, 0, 0,
			A2p1, 0, 0, 0, 0, 0,
			A2s1, 0, 0, 0, 0, 0,
			A31, A32p, A32s, 0, 0, 0,
			A41, A42p, A42s, A43, 0, 0,
			A51, A52p, A52s, A53, A54, 0;

	/* Two photon continuum adds extra spontaneous decay rate to A2s1, but does not produce line radiation */
	_extraAvv << 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			8.23, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0;

	// clang-format on
	// approximation used:
	// A_n,2p = 3 * A_n,2s = 3/4 * (A_n,2 from NIST)

	// The correct way would be to group the levels together in the way described in Hazy II,
	// 6.11. Collision strengths need to be summed, and transition probabilities need to be
	// averaged over. An example grouping would be ( [3d 5/2, 3/2, 1/2] [3p 3/2, 1/2] [3s 1/2] )
	// -> ([2p 3/2, 1/2] [2s 1/2]) Not 100% sure what to to with the multiplicity of the bottom
	// levels though

	DEBUG("Constructed NLevel" << endl);
}

NLevel::NLevel(const Array& frequencyv) : NLevel() { _frequencyv = frequencyv; }

void NLevel::lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const
{
	numLines = 0;
	forAllLinesDo([&](size_t, size_t) { numLines++; });
	lineFreqv.resize(numLines);
	naturalLineWidthv.resize(numLines);
	int index = 0;
	forAllLinesDo([&](size_t upper, size_t lower) {
		lineFreqv[index] = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;
		naturalLineWidthv[index] =
		                (_Avv(upper, lower) + _extraAvv(upper, lower)) / Constant::FPI;
		index++;
	});
}

Eigen::VectorXd NLevel::solveBalance(double atomDensity, double electronDensity, double protonDensity,
                          double temperature, const Array& specificIntensityv, const Array& sourcev,
                          const Array& sinkv)
{
	Eigen::VectorXd nv(_nLv);

	if (specificIntensityv.size() != _frequencyv.size())
		throw range_error("Given ISRF and wavelength vectors do not have the same size");
	if (sourcev.size() != _Ev.size() || sinkv.size() != _Ev.size())
		throw range_error("Source and/or sink term vector(s) of wrong size");

	if (atomDensity > 0)
	{
		Eigen::ArrayXd Cvv(prepareCollisionMatrix(temperature, electronDensity,
		                                          protonDensity));

#ifdef REPORT_LINE_QUALITY
		double maxNorm = 0, minNorm = 1e9;
		for (int lower = 0; lower < _nLv; lower++)
			for (int upper = lower + 1; upper < _nLv; upper++)
				if (_Avv(upper, lower))
				{
					double norm = TemplatedUtils::integrate<double>(
					                _frequencyv, lineProfile(1, 0));
					DEBUG("Line " << upper << " --> " << lower << " has norm "
					              << norm << endl);
					maxNorm = max(norm, maxNorm);
					minNorm = min(norm, minNorm);
				}
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
#endif

		// Calculate BijPij (needs to be redone at each temperature because the line profile
		// can change) Also needs the Cij to calculate collisional broadening
		Eigen::MatrixXd BPvv(prepareAbsorptionMatrix(temperature, specificIntensityv, Cvv));

#ifdef PRINT_MATRICES
		DEBUG("Aij" << endl << _Avv << endl << endl);
		DEBUG("BPij" << endl << BPvv << endl << endl);
		DEBUG("Cij" << endl << Cvv << endl << endl);
#endif
		// Calculate Fij and bi and solve F.n = b
		nv = solveRateEquations(
		                atomDensity, BPvv, Cvv,
		                Eigen::Map<const Eigen::VectorXd>(&sourcev[0], sourcev.size()),
		                Eigen::Map<const Eigen::VectorXd>(&sinkv[0], sinkv.size()), 0);
	}
	return nv;
}

Array NLevel::emissivityv(const Eigen::VectorXd& nv, const Eigen::MatrixXd& Cvv) const
{
	Array total(_frequencyv.size());
	for (int lower = 0; lower < _nLv; lower++)
		for (int upper = lower + 1; upper < _nLv; upper++)
			if (_Avv(upper, lower))
				total += lineIntensityFactor(upper, lower, nv) *
				         lineProfile(upper, lower, nv);
	return total;
}

Array NLevel::opacityv(const Eigen::VectorXd& nv, const Eigen::MatrixXd& Cvv) const
{
	Array total(_frequencyv.size());
	for (int lower = 0; lower < _nLv; lower++)
		for (int upper = lower + 1; upper < _nLv; upper++)
			if (_Avv(upper, lower))
				total += lineOpacityFactor(upper, lower, nv) *
				         lineProfile(upper, lower, nv);
	return total;
}

Array NLevel::scatteringOpacityv(double T, const Eigen::VectorXd& nv, const Eigen::MatrixXd& BPvv,
                                 const Eigen::MatrixXd& Cvv) const
{
	Array total(_frequencyv.size());
	forAllLinesDo([&](size_t upper, size_t lower){
		total += lineOpacityFactor(upper, lower, nv) * lineProfile(upper, lower, T, Cvv) *
		         lineDecayFraction(upper, lower, BPvv, Cvv);};
	return total;
}

void NLevel::forAllLinesDo(function<void(size_t upper, size_t lower)> thingWithLine)
{
	for (int lower = 0; lower < _nLv; lower++)
		for (int upper = lower + 1; upper < _nLv; upper++)
			if (_Avv(upper, lower))
				thingWithLine(upper, lower);
}

void NLevel::forAllLinesDo(function<void(size_t upper, size_t lower)> thingWithLine) const
{
	for (int lower = 0; lower < _nLv; lower++)
		for (int upper = lower + 1; upper < _nLv; upper++)
			if (_Avv(upper, lower))
				thingWithLine(upper, lower);
}

double NLevel::lineIntensityFactor(size_t upper, size_t lower, const Eigen::VectorXd& nv) const
{
	return (_Ev(upper) - _Ev(lower)) / Constant::FPI * nv(upper) * _Avv(upper, lower);
}

double NLevel::lineOpacityFactor(size_t upper, size_t lower, const Eigen::VectorXd& nv) const
{
	double nu_ij = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij /
	                        nu_ij * _Avv(upper, lower);
	double densityFactor = nv(lower) * _gv(upper) / _gv(lower) - nv(upper);
	return constantFactor * densityFactor;
}

Array NLevel::lineProfile(size_t upper, size_t lower, double temperature,
                          const Eigen::MatrixXd& Cvv) const
{
	double nu0 = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;

	double decayRate = _Avv(upper, lower) + _extraAvv(upper, lower) +
	                   Cvv(upper, lower)    // decay rate of top level
	                   + Cvv(lower, upper); // decay rate of bottom level
	// (stimulated emission doesn't count, as it causes no broadening)

	double thermalVelocity = sqrt(Constant::BOLTZMAN * temperature / Constant::HMASS_CGS);

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

double NLevel::lineDecayFraction(size_t upper, size_t lower, const Eigen::MatrixXd& BPvv,
                                 const Eigen::MatrixXd& Cvv) const
{
	return _Avv(upper, lower) / (_Avv.row(upper).sum() + _extraAvv.row(upper).sum() +
	                             BPvv.row(upper).sum() + Cvv.row(upper).sum());
}

Eigen::MatrixXd NLevel::prepareAbsorptionMatrix(double temperature, const Array& specificIntensityv,
                                                const Eigen::MatrixXd& Cvv)
{
	Eigen::MatrixXd BPvv(_nLv, _nLv);
	// Go over the lower triangle, without diagonal
	for (int lower = 0; lower < _nLv; lower++)
	{
		for (int upper = lower + 1; upper < _nLv; upper++)
		{
			// Calculate Pij for the lower triangle (= stimulated emission)
			BPvv(upper, lower) = TemplatedUtils::integrate<double, Array, Array>(
			                _frequencyv, lineProfile(upper, lower, temperature, Cvv) *
			                                             specificIntensityv);

			// Multiply by Bij in terms of Aij, valid for i > j
			double nu_ij = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;
			BPvv(upper, lower) *= Constant::CSQUARE_TWOPLANCK / nu_ij / nu_ij / nu_ij *
			                      _Avv(upper, lower);

			// Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij =
			// Pji
			BPvv(lower, upper) = _gv(upper) / _gv(lower) * BPvv(upper, lower);
		}
		// Set the diagonal to zero
		BPvv(lower, lower) = 0;
	}
	return BPvv;
}

Eigen::MatrixXd NLevel::prepareCollisionMatrix(double temperature, double electronDensity,
                                               double protonDensity)
{
	Eigen::ArrayXd Cvv(_nLv, _nLv);

	auto setElectronCollisionRate = [&](size_t upper, size_t lower, double bigGamma) {
		double kT = Constant::BOLTZMAN * temperature;
		// Equation 6.17 of Hazy II (6.6 Collision strengths)
		Cvv(upper, lower) = bigGamma * 8.6291e-6 / _gv(upper) / sqrt(temperature) *
		                    electronDensity;
		Cvv(lower, upper) = Cvv(upper, lower) * _gv(upper) / _gv(lower) *
		                    exp((_Ev(lower) - _Ev(upper)) / kT);
	};

	// data from 2002-anderson
	// Electron temperatures in electron volt
	vector<double> electronTemperaturesv_eV = {.5, 1., 3., 5., 10., 15., 20., 25.};

	double currentT_eV = Constant::BOLTZMAN * temperature * Constant::ERG_EV;

	// naively inter- and extrapolate this data linearly
	size_t iRight = TemplatedUtils::index(currentT_eV, electronTemperaturesv_eV);
	if (iRight == 0)
		iRight = 1;
	else if (iRight == electronTemperaturesv_eV.size())
		iRight -= 1;
	size_t iLeft = iRight - 1;

	// Effective collision strength
	vector<double> bigGamma2s1v = {2.6e-1,  2.96e-1, 3.26e-1, 3.39e-1,
	                               3.73e-1, 4.06e-1, 4.36e-1, 4.61e-1};
	double bigGamma2s1 = TemplatedUtils::interpolateLinear(
	                currentT_eV, electronTemperaturesv_eV[iLeft],
	                electronTemperaturesv_eV[iRight], bigGamma2s1v[iLeft],
	                bigGamma2s1v[iRight]);
	setElectronCollisionRate(2, 0, bigGamma2s1);

	vector<double> bigGamma2p1v = {4.29e-1, 5.29e-01, 8.53e-01, 1.15e00,
	                               1.81e00, 2.35e00,  2.81e00,  3.20e00};
	double bigGamma2p1 = TemplatedUtils::interpolateLinear(
	                currentT_eV, electronTemperaturesv_eV[iLeft],
	                electronTemperaturesv_eV[iRight], bigGamma2p1v[iLeft],
	                bigGamma2p1v[iRight]);
	setElectronCollisionRate(1, 0, bigGamma2p1);

	// Important for two-photon continuum vs lyman is the l-changing collisions between 2s and
	// 2p 1964-Pengelly eq 43, assuming for qnl = qnl->nl' if l can only be l=1 or l=0
	double mu_m = Constant::HMASS_CGS / 2 / Constant::ELECTRONMASS;
	double constfactor = 9.93e-6 * sqrt(mu_m);

	//(6n^2(n^2 - l^2 - l - 1)) (eq 44)
	double D2p = 24;
	// even though the second term is zero
	double A2p = _Avv.row(index2p).sum() + _extraAvv.row(index2p).sum();
	// (45 should be the correct one for DeltaE <<<)
	double twolog10Rc = 10.95 + log10(temperature / A2p / A2p / mu_m);
	double q2p2s = constfactor * D2p / sqrt(temperature) *
	               (11.54 + log10(temperature / D2p / mu_m) + twolog10Rc);
	Cvv(index2p, index2s) = protonDensity * q2p2s;

	double D2s = 24 * 3;
	double A2s = _Avv.row(index2s).sum() + _extraAvv.row(index2s).sum();
	twolog10Rc = 10.95 + log10(temperature / A2s / A2s / mu_m);
	double q2s2p = constfactor * D2s / sqrt(temperature) *
	               (11.54 + log10(temperature / D2s / mu_m) + twolog10Rc);
	Cvv(index2s, index2p) = protonDensity * q2s2p;

	return Cvv;
}

Eigen::VectorXd NLevel::solveRateEquations(double n, const Eigen::MatrixXd& BPvv,
                                           const Eigen::MatrixXd& Cvv,
                                           const Eigen::VectorXd& sourceTerm,
                                           const Eigen::VectorXd& sinkTerm, int chooseConsvEq)
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	Eigen::MatrixXd Mvv(_Avv.transpose() + _extraAvv.transpose() + BPvv.transpose() +
	                    Cvv.transpose());

	// See equation for Fij (37) in document
	// subtract departure rate from level i to all other levels
	Eigen::MatrixXd departureDiagonal = Mvv.colwise().sum().asDiagonal();
	Mvv -= departureDiagonal;
	Mvv -= sinkTerm.asDiagonal();
	Eigen::VectorXd b(-sourceTerm);

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = Eigen::VectorXd::Ones(Mvv.cols());
	b(chooseConsvEq) = n;

#ifdef PRINT_MATRICES
	DEBUG("System to solve:\n" << Mvv << " * nv\n=\n" << b << endl << endl);
#endif

	// Call the linear solver
	Eigen::VectorXd nv = Mvv.colPivHouseholderQr().solve(b);
	DEBUG("nv" << endl << _nv << endl);

	// Element wise relative errors
	Eigen::ArrayXd diffv = Mvv * nv - b;
	Eigen::ArrayXd errv = diffv / Eigen::ArrayXd(b);
	DEBUG("The relative errors are: " << errv << endl);

	return nv;
}
