#include "NLevel.h"
#include "Constants.h"
#include "SpecialFunctions.h"
#include "TemplatedUtils.h"

using namespace std;

NLevel::NLevel(const Array& frequencyv) : _frequencyv(frequencyv)
{
	// Using NIST data
	// (http://physics.nist.gov/cgi-bin/ASD/lines1.pl?unit=1&line_out=0&bibrefs=1&show_obs_wl=1&show_calc_wl=1&A_out=0&intens_out=1&allowed_out=1&forbid_out=1&conf_out=1&term_out=1&enrg_out=1&J_out=1&g_out=0&spectra=H%20I)
	// 1s, 2p, 2s
	// 3, 4, 5
	_N = 6;
	_Ev << 0., 82258.9191133, 82258.9543992821, 97492.304, 102823.904, 105291.657;
	// Energy in cm^-1, so multiply with hc
	_Ev *= Constant::PLANCKLIGHT;

	_gv << 1, 1, 3, 9, 16, 25;

	_nv = Eigen::VectorXd(_Ev.size());
	_nv(0) = _n;

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
	// A_n,2p = 3 * A_n,2s = 3/4 * A_n,2 from NIST

	_BPvv = Eigen::MatrixXd::Zero(_N);
	_Cvv = Eigen::MatrixXd::Zero(_N);
}

void NLevel::solveBalance(double atomDensity, double electronDensity, double protonDensity,
                          double T, const Array& specificIntensityv, const Array& sourcev,
                          const Array& sinkv)
{
	if (specificIntensityv.size() != _frequencyv.size())
		throw range_error("Given ISRF and wavelength vectors do not have the same size");
	if (sourcev.size() != _Ev.size() || sinkv.size() != _Ev.size())
		throw range_error("Source and/or sink term vector(s) of wrong size");

	_n = atomDensity;
	_ne = electronDensity;
	_np = protonDensity;
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
		prepareAbsorptionMatrix(specificIntensityv);

#ifdef PRINT_MATRICES
		DEBUG("Aij" << endl << _Avv << endl << endl);
		DEBUG("BPij" << endl << _BPvv << endl << endl);
		DEBUG("Cij" << endl << _Cvv << endl << endl);
#endif
		// Calculate Fij and bi and solve F.n = b
		solveRateEquations(Eigen::Map<const Eigen::VectorXd>(&sourcev[0], sourcev.size()),
		                   Eigen::Map<const Eigen::VectorXd>(&sinkv[0], sinkv.size()), 0);
	}
	else
	{
		_nv = Eigen::VectorXd::Zero(_Ev.size());
	}
}

Array NLevel::emissivityv() const
{
	Array total(_frequencyv.size());
	for (size_t lower = 0; lower < _N; lower++)
		for (size_t upper = lower + 1; upper < _N; upper++)
			if (_Avv(upper, lower))
				total += lineIntensityFactor(upper, lower) *
				         lineProfile(upper, lower);
	return total;
}

Array NLevel::opacityv() const
{
	Array total(_frequencyv.size());
	for (size_t lower = 0; lower < _N; lower++)
		for (size_t upper = lower + 1; upper < _N; upper++)
			if (_Avv(upper, lower))
				total += lineOpacityFactor(upper, lower) *
				         lineProfile(upper, lower);
	return total;
}

Array NLevel::scatteringOpacityv() const
{
	Array total(_frequencyv.size());
	for (size_t lower = 0; lower < _N; lower++)
		for (size_t upper = lower + 1; upper < _N; upper++)
			if (_Avv(upper, lower))
				total += lineOpacityFactor(upper, lower) *
				         lineProfile(upper, lower) *
				         lineDecayFraction(upper, lower);
	return total;
}

double NLevel::lineIntensityFactor(size_t upper, size_t lower) const
{
	return (_Ev(upper) - _Ev(lower)) / Constant::FPI * _nv(upper) * _Avv(upper, lower);
}

double NLevel::lineOpacityFactor(size_t upper, size_t lower) const
{
	double nu_ij = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij /
	                        nu_ij * _Avv(upper, lower);
	double densityFactor = _nv(lower) * _gv(upper) / _gv(lower) - _nv(upper);
	return constantFactor * densityFactor;
}

Array NLevel::lineProfile(size_t upper, size_t lower) const
{
	double nu0 = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;

	double decayRate = _Avv(upper, lower) + _extraAvv(upper, lower) +
	                   _Cvv(upper, lower)    // decay rate of top level
	                   + _Cvv(lower, upper); // decay rate of bottom level
	// (stimulated emission doesn't count, as it causes no broadening)

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

double NLevel::lineDecayFraction(size_t upper, size_t lower) const
{
	return _Avv(upper, lower) / (_Avv.row(upper).sum() + _extraAvv.row(upper).sum() +
	                             _BPvv.row(upper).sum() + _Cvv.row(upper).sum());
}

void NLevel::prepareAbsorptionMatrix(const Array& specificIntensityv)
{
	// Go over the lower triangle, without diagonal
	for (size_t lower = 0; lower < _N; lower++)
	{
		for (size_t upper = lower + 1; upper < _N; upper++)
		{
			// Calculate Pij for the lower triangle (= stimulated emission)
			_BPvv(upper, lower) = TemplatedUtils::integrate<double, Array, Array>(
			                _frequencyv,
			                lineProfile(upper, lower) * specificIntensityv);

			// Multiply by Bij in terms of Aij, valid for i > j
			double nu_ij = (_Ev(upper) - _Ev(lower)) / Constant::PLANCK;
			_BPvv(upper, lower) *= Constant::CSQUARE_TWOPLANCK / nu_ij / nu_ij / nu_ij *
			                       _Avv(upper, lower);

			// Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij =
			// Pji
			_BPvv(lower, upper) = _gv(upper) / _gv(lower) * _BPvv(upper, lower);
		}
		// Set the diagonal to zero
		_BPvv(lower, lower) = 0;
	}
}

void NLevel::prepareCollisionMatrix()
{
	// data from 2002-anderson
	vector<double> electronTemperaturesEV = {.5, 1., 3., 5., 10., 15., 20., 25.};

	vector<double> effectiveCollisionStrv2s1 = {2.6e-1,  2.96e-1, 3.26e-1, 3.39e-1,
	                                         3.73e-1, 4.06e-1, 4.36e-1, 4.61e-1};

	vector<double> effectiveCollisionStrv2p1 = {4.29e-1, 5.29e-01, 8.53e-01, 1.15e00,
	                                         1.81e00, 2.35e00,  2.81e00,  3.20e00};
}
