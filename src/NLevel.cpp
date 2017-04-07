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

NLevel::NLevel(const Array& frequencyv, int nLv, const Eigen::VectorXd& ev,
               const Eigen::VectorXd& gv, const Eigen::MatrixXd& avv,
               const Eigen::MatrixXd& extraAvv)
                : _frequencyv(frequencyv), _nLv(nLv), _ev(ev), _gv(gv), _avv(avv),
                  _extraAvv(extraAvv)
{
}

NLevel::NLevel(int nLv, const Eigen::VectorXd& ev, const Eigen::VectorXd& gv,
               const Eigen::MatrixXd& avv, const Eigen::MatrixXd& extraAvv)
                : _nLv(nLv), _ev(ev), _gv(gv), _avv(avv), _extraAvv(extraAvv)
{
}

void NLevel::lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const
{
	numLines = 0;
	forAllLinesDo([&](size_t, size_t) { numLines++; });
	lineFreqv.resize(numLines);
	naturalLineWidthv.resize(numLines);
	int index = 0;
	forAllLinesDo([&](size_t upper, size_t lower) {
		lineFreqv[index] = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
		naturalLineWidthv[index] =
		                (_avv(upper, lower) + _extraAvv(upper, lower)) / Constant::FPI;
		index++;
	});
}

NLevel::Solution NLevel::solveBalance(double atomDensity, double electronDensity,
                                      double protonDensity, double temperature,
                                      const Array& specificIntensityv, const Array& sourcev,
                                      const Array& sinkv) const
{
	Solution s;
	s.n = atomDensity;
	s.T = temperature;
	s.bpvv = Eigen::MatrixXd::Zero(_nLv, _nLv);
	s.cvv = Eigen::MatrixXd::Zero(_nLv, _nLv);
	s.nv = Eigen::VectorXd::Zero(_nLv);

	if (specificIntensityv.size() != _frequencyv.size())
		throw range_error("Given ISRF and wavelength vectors do not have the same size");
	if (sourcev.size() != _ev.size() || sinkv.size() != _ev.size())
		throw range_error("Source and/or sink term vector(s) of wrong size");

	if (atomDensity > 0)
	{
		s.cvv = prepareCollisionMatrix(temperature, electronDensity, protonDensity);
		// Calculate BijPij (needs to be redone at each temperature because the line profile
		// can change) Also needs the Cij to calculate collisional broadening
		s.bpvv = prepareAbsorptionMatrix(specificIntensityv, s.T, s.cvv);

#ifdef REPORT_LINE_QUALITY
		double maxNorm = 0, minNorm = 1e9;
		forAllLinesDo([&](size_t upper, size_t lower) {
			double norm = TemplatedUtils::integrate<double>(
			                _frequencyv, lineProfile(upper, lower, s));
			DEBUG("Line " << upper << " --> " << lower << " has norm " << norm << endl);
			maxNorm = max(norm, maxNorm);
			minNorm = min(norm, minNorm);
		});
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
#endif
#ifdef PRINT_MATRICES
		DEBUG("Aij" << endl << _avv << endl << endl);
		DEBUG("BPij" << endl << s.bpvv << endl << endl);
		DEBUG("Cij" << endl << s.cvv << endl << endl);
#endif
		// Calculate Fij and bi and solve F.n = b
		s.nv = solveRateEquations(
		                s.n, s.bpvv, s.cvv,
		                Eigen::Map<const Eigen::VectorXd>(&sourcev[0], sourcev.size()),
		                Eigen::Map<const Eigen::VectorXd>(&sinkv[0], sinkv.size()), 0);
	}
	return s;
}

Array NLevel::emissivityv(const Solution& info) const
{
	Array total(_frequencyv.size());
	forAllLinesDo([&](size_t upper, size_t lower) {
		total += lineIntensityFactor(upper, lower, info) * lineProfile(upper, lower, info);
	});
	return total;
}

Array NLevel::opacityv(const Solution& info) const
{
	Array total(_frequencyv.size());
	forAllLinesDo([&](size_t upper, size_t lower) {
		total += lineOpacityFactor(upper, lower, info) * lineProfile(upper, lower, info);
	});
	return total;
}

Array NLevel::scatteringOpacityv(const Solution& info) const
{
	Array total(_frequencyv.size());
	forAllLinesDo([&](size_t upper, size_t lower) {
		total += lineOpacityFactor(upper, lower, info) * lineProfile(upper, lower, info) *
		         lineDecayFraction(upper, lower, info);
	});
	return total;
}

void NLevel::forAllLinesDo(function<void(size_t upper, size_t lower)> thingWithLine) const
{
	for (int lower = 0; lower < _nLv; lower++)
		for (int upper = lower + 1; upper < _nLv; upper++)
			if (_avv(upper, lower))
				thingWithLine(upper, lower);
}

double NLevel::lineIntensityFactor(size_t upper, size_t lower, const Solution& info) const
{
	return (_ev(upper) - _ev(lower)) / Constant::FPI * info.nv(upper) * _avv(upper, lower);
}

double NLevel::lineOpacityFactor(size_t upper, size_t lower, const Solution& info) const
{
	double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij /
	                        nu_ij * _avv(upper, lower);
	double densityFactor = info.nv(lower) * _gv(upper) / _gv(lower) - info.nv(upper);
	return constantFactor * densityFactor;
}

inline Array NLevel::lineProfile(size_t upper, size_t lower, const Solution& info) const
{
	return lineProfile(upper, lower, info.T, info.cvv);
}

Array NLevel::lineProfile(size_t upper, size_t lower, double T, const Eigen::MatrixXd& Cvv) const
{
	double nu0 = (_ev(upper) - _ev(lower)) / Constant::PLANCK;

	double decayRate = _avv(upper, lower) + _extraAvv(upper, lower) +
	                   Cvv(upper, lower)    // decay rate of top level
	                   + Cvv(lower, upper); // decay rate of bottom level
	// (stimulated emission doesn't count, as it causes no broadening)

	if (decayRate < 0)
		cout << _avv << endl << _extraAvv << endl << Cvv << endl;

	double thermalVelocity = sqrt(Constant::BOLTZMAN * T / Constant::HMASS_CGS);

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

double NLevel::lineDecayFraction(size_t upper, size_t lower, const Solution& info) const
{
	return _avv(upper, lower) / (_avv.row(upper).sum() + _extraAvv.row(upper).sum() +
	                             info.bpvv.row(upper).sum() + info.cvv.row(upper).sum());
}

Eigen::MatrixXd NLevel::prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
                                                const Eigen::MatrixXd& Cvv) const
{
	Eigen::MatrixXd BPvv = Eigen::MatrixXd::Zero(_nLv, _nLv);
	forAllLinesDo([&](size_t upper, size_t lower) {
		// Calculate Pij for the lower triangle (= stimulated emission)
		BPvv(upper, lower) = TemplatedUtils::integrate<double, Array, Array>(
		                _frequencyv,
		                lineProfile(upper, lower, T, Cvv) * specificIntensityv);

		// Multiply by Bij in terms of Aij, valid for i > j
		double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
		BPvv(upper, lower) *= Constant::CSQUARE_TWOPLANCK / nu_ij / nu_ij / nu_ij *
		                      _avv(upper, lower);

		// Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij =
		// Pji
		BPvv(lower, upper) = _gv(upper) / _gv(lower) * BPvv(upper, lower);
	});
	return BPvv;
}

Eigen::VectorXd NLevel::solveRateEquations(double n, const Eigen::MatrixXd& BPvv,
                                           const Eigen::MatrixXd& Cvv,
                                           const Eigen::VectorXd& sourceTerm,
                                           const Eigen::VectorXd& sinkTerm, int chooseConsvEq) const
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	Eigen::MatrixXd Mvv(_avv.transpose() + _extraAvv.transpose() + BPvv.transpose() +
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
	DEBUG("nv" << endl << nv << endl);

	// Element wise relative errors
	Eigen::ArrayXd diffv = Mvv * nv - b;
	Eigen::ArrayXd errv = diffv / Eigen::ArrayXd(b);
	DEBUG("The relative errors are: " << errv << endl);

	return nv;
}
