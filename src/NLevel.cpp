#include "NLevel.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "LevelDataProvider.h"
#include "TemplatedUtils.h"
#include "Testing.h"

using namespace std;

NLevel::NLevel(shared_ptr<const LevelDataProvider> ldp, const Array& frequencyv, double mass)
		: _ldp(ldp), _frequencyv(frequencyv), _mass(mass), _numLv(_ldp->numLv()),
		  _ev(_ldp->ev()), _gv(_ldp->gv()), _avv(_ldp->avv()),
		  _extraAvv(_ldp->extraAvv())
{
	// Do a sanity check: All active transitions must be downward ones in energy
	forActiveLinesDo([&](size_t upper, size_t lower) {
		if (_ev(upper) - _ev(lower) < 0)
		{
			cout << "Eu - El = " << _ev(upper) - _ev(lower) << endl;
			cout << "u l = " << upper << " " << lower
			     << " Aul = " << _avv(upper, lower) << endl;
			Error::runtime("There is an upward A-coefficient. This can't be "
				       "correct");
		}
	});
}

NLevel::~NLevel() = default;

void NLevel::lineInfo(int& numLines, Array& lineFreqv, Array& naturalLineWidthv) const
{
	// Count the number of active lines
	numLines = 0;
	forActiveLinesDo([&](size_t, size_t) { numLines++; });
	lineFreqv.resize(numLines);
	naturalLineWidthv.resize(numLines);

	// Calculate the natural line width for these transitions
	int index = 0;
	forActiveLinesDo([&](size_t upper, size_t lower) {
		lineFreqv[index] = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
		naturalLineWidthv[index] =
				(_avv(upper, lower) + _extraAvv(upper, lower)) / Constant::FPI;
		index++;
	});
}

NLevel::Solution NLevel::solveBalance(double density, const EVector& speciesNv,
				      double temperature, const Spectrum& specificIntensity,
				      const EVector& sourcev, const EVector& sinkv) const
{
	Solution s;
	s.n = density;
	s.T = temperature;
	s.bpvv = EMatrix::Zero(_numLv, _numLv);
	s.cvv = EMatrix::Zero(_numLv, _numLv);
	s.nv = EVector::Zero(_numLv);

	if (density > 0)
	{
		s.cvv = _ldp->cvv(temperature, speciesNv);
		/* Calculate BijPij (needs to be redone at each temperature because the line
		   profile can change). Also needs the Cij to calculate collisional
		   broadening. */
		s.bpvv = prepareAbsorptionMatrix(specificIntensity, s.T, s.cvv);
// #define REPORT_LINE_QUALITY
#ifdef REPORT_LINE_QUALITY
		// Full integral of the line profile, to check the discretization of the output
		// (emission) grid.
		double maxNorm = 0, minNorm = 1e9;
		forActiveLinesDo([&](size_t upper, size_t lower) {
			double norm = TemplatedUtils::integrate<double>(
					_frequencyv, lineProfile_array(upper, lower, s));
			DEBUG("lineProfile_array " << upper << " --> " << lower << " has norm "
						   << norm << endl);
			maxNorm = max(norm, maxNorm);
			minNorm = min(norm, minNorm);
		});
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
		maxNorm = 0;
		minNorm = 1e9; // any large number will do
		// Integral over contant spectrum, using the LineProfile class. An integration
		// grid is chosen internally, and we check the quality of it here (at least for
		// now.)
		double minFreq = _frequencyv[0];
		double maxFreq = _frequencyv[_frequencyv.size() - 1];
		Array someFreqs = {minFreq, (minFreq + maxFreq) / 2, maxFreq};
		Spectrum flat(someFreqs, Array(1, someFreqs.size()));
		forActiveLinesDo([&](size_t upper, size_t lower) {
			auto lp = lineProfile(upper, lower, s);
			double norm = lp.integrateSpectrum(flat);
			DEBUG("LineProfile " << upper << " --> " << lower << " has norm "
					     << norm << endl);
			maxNorm = max(norm, maxNorm);
			minNorm = min(norm, minNorm);
		});
		DEBUG("Max profile norm = " << maxNorm << endl);
		DEBUG("Min profile norm = " << minNorm << endl);
#endif
#ifdef PRINT_LEVEL_MATRICES
		DEBUG("Aij" << endl << _avv << endl << endl);
		DEBUG("BPij" << endl << s.bpvv << endl << endl);
		DEBUG("Cij" << endl << s.cvv << endl << endl);
#endif
		// Calculate Fij and bi and solve F.n = b
		s.nv = solveRateEquations(s.n, s.bpvv, s.cvv, sourcev, sinkv, 0);
	}
	return s;
}

NLevel::Solution NLevel::solveLTE(double density, const EVector& speciesNv, double T,
				  const Spectrum& specificIntensity) const
{
	NLevel::Solution s;
	s.n = density;
	s.T = T;
	s.nv = density * solveBoltzmanEquations(T);
	s.cvv = _ldp->cvv(T, speciesNv);
	s.bpvv = prepareAbsorptionMatrix(specificIntensity, T, s.cvv);
	return s;
}

NLevel::Solution NLevel::solveZero(double T) const
{
	NLevel::Solution s;
	s.n = 0;
	s.T = T;
	s.nv = EVector::Zero(_numLv);
	s.cvv = EMatrix::Zero(_numLv, _numLv);
	s.bpvv = EMatrix::Zero(_numLv, _numLv);
	return s;
}

Array NLevel::emissivityv(const Solution& s) const { return lineEmissivityv(s); }

Array NLevel::lineEmissivityv(const Solution& s) const
{
	Array total(_frequencyv.size());
	forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = lineIntensityFactor(upper, lower, s);
		LineProfile lp = lineProfile(upper, lower, s);
		lp.addToSpectrum(_frequencyv, total, factor);
	});
	return total;
}

Array NLevel::opacityv(const Solution& s) const { return lineOpacityv(s); }

Array NLevel::lineOpacityv(const Solution& s) const
{
	Array total(_frequencyv.size());
	forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = lineOpacityFactor(upper, lower, s);
		LineProfile lp = lineProfile(upper, lower, s);
		lp.addToSpectrum(_frequencyv, total, factor);
	});
	return total;
}

double NLevel::heating(const Solution& s) const
{
	double total = 0;
	for (size_t ini = 0; ini < _numLv; ini++)
	{
		for (size_t fin = 0; fin < _numLv; fin++)
		{
			// Downward transitions inject kinetic energy into the medium
			if (_ev(ini) > _ev(fin))
			{
				double cul = s.cvv(ini, fin);
				if (cul > 0)
					total += (_ev(ini) - _ev(fin)) * cul * s.nv(ini);
			}
		}
	}
	return total;
}

double NLevel::cooling(const Solution& s) const
{
	double total = 0;
	for (size_t ini = 0; ini < _numLv; ini++)
	{
		for (size_t fin = 0; fin < _numLv; fin++)
		{
			// Upward transitions absorb kinetic energy from the medium
			if (_ev(ini) < _ev(fin))
			{
				double clu = s.cvv(ini, fin);
				if (clu > 0)
					total += (_ev(fin) - _ev(ini)) * clu * s.nv(ini);
			}
		}
	}
	return total;
}

EMatrix NLevel::prepareAbsorptionMatrix(const Spectrum& specificIntensity, double T,
					const EMatrix& Cvv) const
{
	EMatrix BPvv = EMatrix::Zero(_numLv, _numLv);
	const Array& v = specificIntensity.valuev();
	auto maxIt = max_element(begin(v), end(v));
	double spectrumMax = *maxIt;

	Array highResIv(_frequencyv.size());
	for (size_t i = 0; i < _frequencyv.size(); i++)
		highResIv[i] = specificIntensity.evaluate(_frequencyv[i]);
	Spectrum highResSpec(_frequencyv, highResIv);

	forActiveLinesDo([&](size_t upper, size_t lower) {
		// Calculate Pij for the lower triangle (= stimulated emission)
		LineProfile lp = lineProfile(upper, lower, T, Cvv);
		double lowResIntegral = lp.integrateSpectrum(specificIntensity);
// #define REPORT_SPEC_INTEGRAL
#ifdef REPORT_SPEC_INTEGRAL
		double highResIntegral = lp.integrateSpectrum(highResSpec);

		const Array& lpav = lineProfile_array(upper, lower, T, Cvv);
		double manualIntegral = TemplatedUtils::integrate<double, Array, Array>(
				_frequencyv, lpav * highResIv);

		double hrRatio = highResIntegral / lowResIntegral;
		double mnRatio = manualIntegral / lowResIntegral;

		if (abs(hrRatio - 1.) > 1.e-6 || abs(mnRatio - 1.) > 1.e-6)
			cout << lowResIntegral << "\t HR: " << hrRatio << "\t MR:" << mnRatio
			     << endl;
#endif /* REPORT_SPEC_INTEGRAL */
		BPvv(upper, lower) = lowResIntegral;

		// Multiply by Bij in terms of Aij, valid for i > j
		double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
		BPvv(upper, lower) *= Constant::CSQUARE_TWOPLANCK / nu_ij / nu_ij / nu_ij *
				      _avv(upper, lower);

		/* Derive the upper triangle (= absorption) using gi Bij = gj Bji and Pij =
		   Pji. */
		BPvv(lower, upper) = _gv(upper) / _gv(lower) * BPvv(upper, lower);
	});
	return BPvv;
}

EMatrix NLevel::netTransitionRate(const EMatrix& BPvv, const EMatrix& Cvv) const
{
	return _avv.transpose() + _extraAvv.transpose() + BPvv.transpose() + Cvv.transpose();
}

EVector NLevel::solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
				   const EVector& sourcev, const EVector& sinkv,
				   int chooseConsvEq) const
{
	// Initialize Mij as Aji + PBji + Cji
	// = arrival rate in level i from level j
	EMatrix Mvv = netTransitionRate(BPvv, Cvv);

	// See equation for Fij in gasphysics document
	// Subtract departure rate from level i to all other levels
	EMatrix departureDiagonal = Mvv.colwise().sum().asDiagonal();
	Mvv -= departureDiagonal; // Fij
	Mvv -= sinkv.asDiagonal(); // Fij - diag[d_i]
	EVector f(-sourcev); // -f_i

	// Replace row by a conservation equation
	Mvv.row(chooseConsvEq) = EVector::Ones(Mvv.cols());
	f(chooseConsvEq) = n;

#ifdef PRINT_LEVEL_MATRICES
	DEBUG("System to solve:\n" << Mvv << " * nv\n=\n" << f << endl << endl);
#endif
	// Call the linear solver for the system sum_j (Fij - diag[d_i]) * n_j = -f_i
	EVector nv = Mvv.colPivHouseholderQr().solve(f);
	// DEBUG("nv" << endl << nv << endl);

	// Hack: put populations = 0 if they were negative due to precision issues
	nv = nv.array().max(0);

	// Element wise relative errors
	EArray diffv = Mvv * nv - f;
	EArray errv = diffv / EArray(f);
	// DEBUG("The relative errors are:\n" << errv << endl);

	return nv;
}

EVector NLevel::solveBoltzmanEquations(double T) const
{
	// Degeneracy factor
	EVector pv{_gv};
	double pSum{0};
	double kT = Constant::BOLTZMAN * T;
	for (int i = 0; i < _ev.size(); i++)
	{
		// Exponential factor
		pv(i) *= exp(-_ev(i) / kT);
		// Normalization
		pSum += pv(i);
	}
	return pv / pSum;
}

void NLevel::forActiveLinesDo(function<void(size_t ini, size_t fin)> thing) const
{
	// Execute the same function for all transitions that are optically active.
	for (size_t fin = 0; fin < _numLv; fin++)
		for (size_t ini = 0; ini < _numLv; ini++)
			if (_avv(ini, fin))
				thing(ini, fin);
}

double NLevel::lineIntensityFactor(size_t upper, size_t lower, const Solution& s) const
{
	return (_ev(upper) - _ev(lower)) / Constant::FPI * s.nv(upper) * _avv(upper, lower);
}

double NLevel::lineOpacityFactor(size_t upper, size_t lower, const Solution& s) const
{
	double nu_ij = (_ev(upper) - _ev(lower)) / Constant::PLANCK;
	double constantFactor = Constant::LIGHT * Constant::LIGHT / 8. / Constant::PI / nu_ij /
				nu_ij * _avv(upper, lower);
	double densityFactor = s.nv(lower) * _gv(upper) / _gv(lower) - s.nv(upper);
	double result = constantFactor * densityFactor;
#ifdef SANITY
	if (result < 0)
		result = 0;
#endif
	return result;
}

LineProfile NLevel::lineProfile(size_t upper, size_t lower, const Solution& s) const
{
	return lineProfile(upper, lower, s.T, s.cvv);
}

LineProfile NLevel::lineProfile(size_t upper, size_t lower, double T, const EMatrix& Cvv) const
{
	double nu0 = (_ev(upper) - _ev(lower)) / Constant::PLANCK;

	double decayRate = _avv(upper, lower) + _extraAvv(upper, lower) +
			   Cvv(upper, lower) // decay rate of top level
			   + Cvv(lower, upper); // decay rate of bottom level
	// (stimulated emission doesn't count, as it causes no broadening)

	double thermalVelocity = sqrt(Constant::BOLTZMAN * T / _mass);

	// Half the FWHM of the Lorentz
	double halfWidth = decayRate / Constant::FPI;

	// The standard deviation in frequency units. It is about half of the FWHM for a
	// Gaussian
	double sigma_nu = nu0 * thermalVelocity / Constant::LIGHT;

	return LineProfile(nu0, sigma_nu, halfWidth);
}

Array NLevel::lineProfile_array(size_t upper, size_t lower, const Solution& s) const
{
	return lineProfile_array(upper, lower, s.T, s.cvv);
}

Array NLevel::lineProfile_array(size_t upper, size_t lower, double T, const EMatrix& Cvv) const
{
	Array result(_frequencyv.size());
	LineProfile lp = lineProfile(upper, lower, T, Cvv);
	for (size_t n = 0; n < _frequencyv.size(); n++)
		result[n] = lp(_frequencyv[n]);
	return result;
}
