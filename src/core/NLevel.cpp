#include "NLevel.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GasStruct.hpp"
#include "LevelDataProvider.hpp"
#include "Options.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"

using namespace std;

NLevel::NLevel(shared_ptr<const LevelDataProvider> ldp, double mass)
                : _ldp(ldp), _mass(mass), _numLv(_ldp->numLv()), _ev(_ldp->ev()),
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

EMatrix NLevel::totalTransitionRatesvv(const Spectrum& specificIntensity, const GasStruct& gas,
                                       EMatrix* cvv_p) const
{
	EMatrix cvv = _ldp->cvv(gas);
	if (cvv_p)
		*cvv_p = cvv;
	EMatrix bpvv = prepareAbsorptionMatrix(specificIntensity, gas._T, cvv);
	if (Options::nlevel_printLevelMatrices)
	{
		DEBUG("Aij" << endl << _avv << endl << endl);
		DEBUG("BPij" << endl << bpvv << endl << endl);
		DEBUG("Cij" << endl << cvv << endl << endl);
	}
	return _avv + _extraAvv + bpvv + cvv;
}

NLevel::Solution NLevel::solveLTE(double density, const GasStruct& gas) const
{
	NLevel::Solution s;
	s.n = density;
	s.T = gas._T;
	s.nv = density * solveBoltzmanEquations(gas._T);
	s.cvv = _ldp->cvv(gas);
	return s;
}

NLevel::Solution NLevel::solveZero(double T) const
{
	NLevel::Solution s;
	s.n = 0;
	s.T = T;
	s.nv = EVector::Zero(_numLv);
	s.cvv = EMatrix::Zero(_numLv, _numLv);
	return s;
}

Array NLevel::emissivityv(const Solution& s, const Array& eFrequencyv) const
{
	return lineEmissivityv(s, eFrequencyv);
}

Array NLevel::lineEmissivityv(const Solution& s, const Array& eFrequencyv) const
{
	Array total(eFrequencyv.size());
	forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = lineIntensityFactor(upper, lower, s);
		LineProfile lp = lineProfile(upper, lower, s);
		lp.addToBinned(eFrequencyv, total, factor);
	});
	return total;
}

Array NLevel::opacityv(const Solution& s, const Array& oFrequencyv) const
{
	return lineOpacityv(s, oFrequencyv);
}

Array NLevel::lineOpacityv(const Solution& s, const Array& oFrequencyv) const
{
	Array total(oFrequencyv.size());
	forActiveLinesDo([&](size_t upper, size_t lower) {
		double factor = lineOpacityFactor(upper, lower, s);
		LineProfile lp = lineProfile(upper, lower, s);
		// lp.addToSpectrum(oFrequencyv, total, factor);
		lp.addToBinned(oFrequencyv, total, factor);
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

double NLevel::netheating(const Solution& s) const
{
	double total = 0;
	for (size_t up = 0; up < _numLv; up++)
	{
		for (size_t lo = 0; lo < up; lo++)
		{
			double cul = s.cvv(up, lo);
			double clu = s.cvv(lo, up);
			if (cul > 0)
				total += (_ev(up) - _ev(lo)) *
				         (cul * s.nv(up) - clu * s.nv(lo));
		}
	}
	return total;
}

EMatrix NLevel::prepareAbsorptionMatrix(const Spectrum& specificIntensity, double T,
                                        const EMatrix& Cvv) const
{
	EMatrix BPvv = EMatrix::Zero(_numLv, _numLv);
	double spectrumMax = specificIntensity.valMax();
	forActiveLinesDo([&](size_t upper, size_t lower) {
		// Calculate Pij for the lower triangle (= stimulated emission)
		LineProfile lp = lineProfile(upper, lower, T, Cvv);
		double lowResIntegral = lp.integrateSpectrum(specificIntensity, spectrumMax);
// #define REPORT_SPEC_INTEGRAL
#ifdef REPORT_SPEC_INTEGRAL
		auto f = [&](double x) -> double {
			return lp(x) * specificIntensity.evaluate(x);
		};
		// TODO: better use log grid here, and put this in a separate test function
		size_t many_points = 1e6;
		double manualIntegral = TemplatedUtils::integrateFunction<double>(
		                f, specificIntensity.freqMin(), specificIntensity.freqMax(),
		                many_points);

		double ratio = manualIntegral / lowResIntegral;

		if (abs(ratio - 1.) > 1.e-6)
			cout << lowResIntegral << "\t MR:" << ratio << endl;
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

EVector NLevel::solveBoltzmanEquations(double T) const
{
	double eMin = _ev.minCoeff();
	// Degeneracy factor
	EVector pv{_gv};
	// Partition function
	double pSum{0};
	double kT = Constant::BOLTZMAN * T;
	for (int i = 0; i < _ev.size(); i++)
	{
		// Exponential factor
		pv(i) *= exp((eMin - _ev(i)) / kT);
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
