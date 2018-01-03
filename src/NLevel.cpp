#include "NLevel.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "LevelDataProvider.h"
#include "TemplatedUtils.h"
#include "Testing.h"

using namespace std;

NLevel::NLevel(shared_ptr<const LevelDataProvider> ldp, const Array& frequencyv)
                : _ldp(ldp), _frequencyv(frequencyv), _numLv(_ldp->numLv()), _ev(_ldp->ev()),
                  _gv(_ldp->gv()), _avv(_ldp->avv()), _extraAvv(_ldp->extraAvv())
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
                                      double temperature, const Array& specificIntensityv,
                                      const EVector& sourcev, const EVector& sinkv) const
{
	Solution s;
	s.n = density;
	s.T = temperature;
	s.bpvv = EMatrix::Zero(_numLv, _numLv);
	s.cvv = EMatrix::Zero(_numLv, _numLv);
	s.nv = EVector::Zero(_numLv);

	if (specificIntensityv.size() != _frequencyv.size())
		Error::runtime("Given ISRF and frequency vectors do not have the same size " +
		               std::to_string(specificIntensityv.size()) + " vs " +
		               std::to_string(_frequencyv.size()));

	if (density > 0)
	{
		s.cvv = _ldp->cvv(temperature, speciesNv);
		/* Calculate BijPij (needs to be redone at each temperature because the line
		   profile can change) Also needs the Cij to calculate collisional broadening. */
		s.bpvv = prepareAbsorptionMatrix(specificIntensityv, s.T, s.cvv);

#ifdef REPORT_LINE_QUALITY
		double maxNorm = 0, minNorm = 1e9;
		forActiveLinesDo([&](size_t upper, size_t lower) {
			double norm = TemplatedUtils::integrate<double>(
			                _frequencyv, lineProfile(upper, lower, s));
			DEBUG("Line " << upper << " --> " << lower << " has norm " << norm
			              << endl);
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
                                  const Array& specificIntensityv) const
{
	NLevel::Solution s;
	s.n = density;
	s.T = T;
	s.nv = density * solveBoltzmanEquations(T);
	s.cvv = _ldp->cvv(T, speciesNv);
	s.bpvv = prepareAbsorptionMatrix(specificIntensityv, T, s.cvv);
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
		addLine(total, upper, lower, factor, s.T, s.cvv);
	});
	return total;
}

Array NLevel::opacityv(const Solution& s) const { return lineOpacityv(s); }

Array NLevel::lineOpacityv(const Solution& s) const
{
	Array total(_frequencyv.size());
	forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = lineOpacityFactor(upper, lower, s);
		addLine(total, upper, lower, factor, s.T, s.cvv);
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

EMatrix NLevel::prepareAbsorptionMatrix(const Array& specificIntensityv, double T,
                                        const EMatrix& Cvv) const
{
	EMatrix BPvv = EMatrix::Zero(_numLv, _numLv);
	forActiveLinesDo([&](size_t upper, size_t lower) {
		// Calculate Pij for the lower triangle (= stimulated emission)
		BPvv(upper, lower) = TemplatedUtils::integrate<double, Array, Array>(
		                _frequencyv,
		                lineProfile_array(upper, lower, T, Cvv) * specificIntensityv);

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

	// See equation for Fij (37) in document
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
	DEBUG("nv" << endl << nv << endl);

	// Hack: put populations = 0 if they were negative due to precision issues
	nv = nv.array().max(0);

	// Element wise relative errors
	EArray diffv = Mvv * nv - f;
	EArray errv = diffv / EArray(f);
	DEBUG("The relative errors are:\n" << errv << endl);

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

	double thermalVelocity = sqrt(Constant::BOLTZMAN * T / Constant::HMASS_CGS);

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

void NLevel::addLine(Array& spectrumv, size_t upper, size_t lower, double factor, double T,
                     const EMatrix& Cvv) const
{
	LineProfile lp = lineProfile(upper, lower, T, Cvv);

#define OPTIMIZED_LINE_SPECTRUM
#ifdef OPTIMIZED_LINE_SPECTRUM
	// Only calculate the voigt function for significant contributions of this line to the
	// total spectrum. Start at the line center (assumes there is a suitable frequency point
	// in the grid)
	size_t iCenter = TemplatedUtils::index(nu0, _frequencyv);
	const double CUTOFFWINGCONTRIBUTION = 1e-9;
	double wingThres = abs(1e-6 * lineValuef(iCenter));

	// Add values for center and right wing:
	for (size_t i = iCenter; i < _frequencyv.size(); i++)
	{
		double value = factor * lp(_frequencyv[i]);
		spectrumv[i] += value;

		// Once we are in the wing (arbitrarily defined here), check if we can start
		// ignoring some points. Stop evaluating the wing if its contribution becomes
		// negligible
		double absval = abs(value);
		if (absval < wingThres && absval < abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
		{
			// Keep skipping points until the spectrum drops below the threshold for
			// significance of the line again, or until we reach the end of the
			// spectrum.
			while (true)
			{
				// If we reach te end of the spectrum, break. The main loop will
				// exit.
				if (i >= _frequencyv.size())
					break;
				// If the current value is significant compared to the spectrum
				// at index i, roll back by 1 and break. The main loop will then
				// increment back to the current i, and calculate the value of
				// the line profile for it.
				if (absval > abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
				{
					i--;
					break;
				}
				// If the current value is still too small compared to the
				// current spectrum at i (remember that our wing is descending
				// and becoming even smaller), skip this point.
				i++;
			}
		}
	}

	// Add values for left wing
	if (iCenter > 0)
	{
		// Stops when i == 0 (i-- decrements to i - 1, but returns the original, so the
		// loop safely stops when the body exits with i == 0)
		for (size_t i = iCenter; i-- > 0;)
		{
			double value = factor * lp(_frequencyv[i]);
			spectrumv[i] += value;
			double absval = abs(value);
			// For an insignificant wing value
			if (value < wingThres &&
			    absval < abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
			{
				// Move through the spectrum until it the value of the spectrum
				// is small enough for the value of the wing point to be
				// significant again.
				while (true)
				{
					if (i == 0)
						break;
					if (absval > abs(spectrumv[i] * CUTOFFWINGCONTRIBUTION))
					{
						// Found a significant value! Go back up, then
						// let the main loop continue.
						i++;
						break;
					}
					// Did not find significant value, move further in the
					// left wing.
					i--;
				}
			}
		}
	}
#else
	// Add the whole line to the spectrum
	for (size_t i = 0; i < _frequencyv.size(); i++)
		spectrumv[i] += factor * lp(_frequencyv[i]);
#endif
}
