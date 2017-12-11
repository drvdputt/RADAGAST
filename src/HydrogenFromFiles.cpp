#include "HydrogenFromFiles.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "Error.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"

using namespace std;

namespace
{
const int cNMAX = 5;

inline vector<int> twoJplus1range(int l)
{
	/* If l > 0, then j can be either l + 0.5 or l - 0.5. If l==0, j is always 1/2 and
	   2j+1 can only be 2 */
	return l > 0 ? vector<int>({2 * l, 2 * l + 2}) : vector<int>({2});
}
} /* namespace */

HydrogenFromFiles::HydrogenFromFiles(int resolvedUpTo)
                : _resolvedUpTo(resolvedUpTo), _ine{SpeciesIndex::ine()},
                  _inp{SpeciesIndex::inp()}
{
	if (_resolvedUpTo > cNMAX)
		Error::rangeCheck<int>("Number of resolved levels", _resolvedUpTo, 2, cNMAX);
	// Read-in, and steps that are safe during read-in
	readData();
	// Steps that need to happen after read-in
	prepareForOutput();
}

void HydrogenFromFiles::readData()
{
	const string basename{REPOROOT "/dat/CHIANTI_8.0.6_data/h/h_1/h_1"};

	//-----------------//
	// READ LEVEL DATA //
	//-----------------//
	ifstream elvlc = IOTools::ifstreamFile(basename + ".elvlc");
	string line;

	// Start with the first line
	getline(elvlc, line);
	while (line.compare(1, 2, "-1"))
	{
		// Read the different parts of the line
		int lvIndex, twoSplus1;
		string config;
		char lSymbol;
		double j, observedEnergy, theoreticalEnergy;
		istringstream(line) >> lvIndex >> config >> twoSplus1 >> lSymbol >> j >>
		                observedEnergy >> theoreticalEnergy;
#ifdef ECHO_READIN
		DEBUG(lvIndex << " " << config << " " << twoSplus1 << " " << lSymbol << " " << j
		              << " " << observedEnergy << " " << theoreticalEnergy << " "
		              << endl);
#endif

		// Get the first number from the config string
		int n;
		istringstream(config) >> n;

		// Translate the angular momentum letter
		int l = _lNumberm.at(lSymbol);

		// Store 2j+1 (static cast to make it clear that j is not integer)
		int twoJplus1 = static_cast<int>(2 * j + 1);

		// Convert the energy from cm-1 to erg
		double e = observedEnergy * Constant::LIGHT * Constant::PLANCK;
#ifdef ECHO_READIN
		DEBUG("n " << n << " l " << l << " 2j+1 " << twoJplus1 << " e "
		           << e * Constant::ERG_EV << endl);
#endif

		/* The level indices in the data structures will go from 0 to number of levels minus
		   one. The quantum numbers are also used as keys in a map, so we can quickly
		   retrieve the index for a given configuration.  The level indices in the file go
		   from 1 to the number of levels, while those in the map and in all the vectors
		   will go from 0 to numL - 1. */
		_chiantiLevelv.emplace_back(n, l, twoJplus1, e);
		_nljToChiantiIndexm.insert({{n, l, twoJplus1}, lvIndex - 1});

		// Go the the next line
		getline(elvlc, line);
	}
	elvlc.close();
	_chiantiNumLvl = _chiantiLevelv.size();

	//-----------------//
	// READ EINSTEIN A //
	//-----------------//
	_chiantiAvv = EMatrix::Zero(_chiantiNumLvl, _chiantiNumLvl);
	ifstream wgfa = IOTools::ifstreamFile(basename + ".wgfa");
	getline(wgfa, line);
	while (line.compare(1, 2, "-1"))
	{
		int leftIndex, rightIndex;
		double wavAngstrom, gf, A;
		istringstream(line) >> leftIndex >> rightIndex >> wavAngstrom >> gf >> A;

		/* A comment in the cloudy code recommended to do this, as there are apparently some
		   files in the CHIANTI database where the left index represents the upper level of
		   the transition: */
		int upperIndex = max(leftIndex, rightIndex);
		int lowerIndex = min(leftIndex, rightIndex);
#ifdef ECHO_READIN
		DEBUG(lowerIndex << " " << upperIndex << " " << wavAngstrom << " " << gf << " "
		                 << A << endl);
#endif
		// Zero means two-photon transition, see CHIANTI user guide.
		_chiantiAvv(upperIndex - 1, lowerIndex - 1) = wavAngstrom > 0 ? A : 0;

		getline(wgfa, line);
	}
	wgfa.close();
#ifdef ECHO_READIN
	DEBUG("All chianti A coefficients:" << endl);
	DEBUG(_chiantiAvv << endl);
#endif
	//---------------------//
	// READ COLLISION DATA //
	//---------------------//
	int andersonIndex = 1;
	for (int n = 0; n <= cNMAX; n++)
	{
		for (int l = 0; l < n; l++)
		{
			_nlToAndersonIndexm.insert({{n, l}, andersonIndex});
			andersonIndex++;
		}
	}

	ifstream h_coll_str = IOTools::ifstreamRepoFile("dat/h_coll_str.dat");
	getline(h_coll_str, line);
	getline(h_coll_str, line);
	while (line.compare(0, 2, "-1"))
	{
		istringstream iss = istringstream(line);
		int upperIndex, lowerIndex;
		iss >> upperIndex >> lowerIndex;
		Array Upsilonv(8);
		for (int t = 0; t < 8; t++)
		{
			iss >> Upsilonv[t];
#ifdef ECHO_READIN
			DEBUG("temp" << _andersonTempv[t]);
			DEBUG(" " << Upsilonv[t] << " ");
#endif
		}
#ifdef ECHO_READIN
		DEBUG(endl);
#endif
		_andersonUpsilonvm[{upperIndex, lowerIndex}] = Upsilonv;
		getline(h_coll_str, line);
	}
	h_coll_str.close();
}

void HydrogenFromFiles::prepareForOutput()
{
	//-------------------------------------//
	// SET PARAMETERS FOR HELP WITH OUTPUT //
	//-------------------------------------//
	_levelOrdering.clear();
	int n = 1;
	while (n <= _resolvedUpTo)
	{
		for (int l = 0; l < n; l++)
		{
			_levelOrdering.emplace_back(n, l, energy(n, l));
			_nlToOutputIndexm.insert({{n, l}, _levelOrdering.size() - 1});
		}
		n++;
	}
	while (n <= cNMAX)
	{
		_levelOrdering.emplace_back(n, energy(n));
		_nlToOutputIndexm.insert({{n, -1}, _levelOrdering.size() - 1});
		n++;
	}
	_numL = _levelOrdering.size();

	_totalAv = EVector::Zero(_numL);
	for (size_t i = 0; i < _numL; i++)
		for (size_t f = 0; f < _numL; f++)
			_totalAv[i] += einsteinA(_levelOrdering[i], _levelOrdering[f]);
}

size_t HydrogenFromFiles::numLv() const { return _numL; }

size_t HydrogenFromFiles::indexOutput(int n, int l) const
{
	return _nlToOutputIndexm.at({n, l});
}

EVector HydrogenFromFiles::ev() const
{
	EVector the_ev(_numL);
	for (size_t i = 0; i < _numL; i++)
		the_ev(i) = _levelOrdering[i].e();
	return the_ev;
}

EVector HydrogenFromFiles::gv() const
{
	EVector the_gv(_numL);
	for (size_t i = 0; i < _numL; i++)
		the_gv(i) = _levelOrdering[i].g();
	return the_gv;
}

EMatrix HydrogenFromFiles::avv() const
{
	EMatrix the_avv(_numL, _numL);
	for (size_t i = 0; i < _numL; i++)
	{
		const HydrogenLevel& initial = _levelOrdering[i];
		int ni = initial.n();
		int li = initial.l();
		for (size_t f = 0; f < _numL; f++)
		{
			const HydrogenLevel& final = _levelOrdering[f];
			int nf = final.n();
			int lf = final.l();
			// Both resolved
			if (initial.nlResolved() && final.nlResolved())
				the_avv(i, f) = einsteinA(ni, li, nf, lf);
			// Collapsed to resolved
			else if (initial.nCollapsed() && final.nlResolved())
				the_avv(i, f) = einsteinA(ni, nf, lf);
			// Resolved to collapsed should not exist (downward) or be zero (upward).
			else if (initial.nlResolved() && final.nCollapsed())
			{
				/* The collapsed-collapsed equivalent should always be 0 in this
				   case. */
				the_avv(i, f) = einsteinA(ni, nf);
				assert(the_avv(i, f) == 0.);
			}
			// Collapsed to collapsed
			else
				the_avv(i, f) = einsteinA(ni, nf);
		}
	}
#ifdef PRINT_LEVEL_MATRICES
	DEBUG("Einstein A matrix:" << endl);
	DEBUG(the_avv << endl);
#endif
	return the_avv;
}

EMatrix HydrogenFromFiles::extraAvv() const
{
	EMatrix the_extra = EMatrix::Zero(_numL, _numL);

	size_t index1s = indexOutput(1, 0);

	// Check if the n2 level is resolved, and retrieve the Key, Value pair
	auto index2sIt = _nlToOutputIndexm.find({2, 0});
	/* If the level is collapsed, the pair 0,2 won't be found, and the find function will return
	   'end'. So if the result is not 'end', this means that the resolved 2,0 level was
	   found. */
	bool n2Resolved = index2sIt != _nlToOutputIndexm.end();

	// Hardcode the two-photon decay of 2s to 1s
	double rate = 8.229;
	// If 2nd level is resolved, just put the value in the right place
	if (n2Resolved)
	{
		size_t index2s = index2sIt->second;
		the_extra(index2s, index1s) = rate;
	}
	/* If it is collapsed, assume it is well mixed. So we take the average of the two-photon
	   rates from 2s (8.229) and 2p (0). */
	else
	{
		size_t index2 = indexOutput(2, -1);
		the_extra(index2, index1s) = rate / 4.; // + 0 * 3 / 4.
	}
#ifdef PRINT_LEVEL_MATRICES
	DEBUG("Extra A" << endl);
	DEBUG(the_extra << endl);
#endif
	return the_extra;
}

array<size_t, 2> HydrogenFromFiles::twoPhotonIndices() const
{
	// If any of the levels is not resolved on l, just return the index of the collapsed level.
	size_t upper = _resolvedUpTo >= 2 ? indexOutput(2, 0) : indexOutput(2, -1);
	size_t lower = _resolvedUpTo >= 1 ? indexOutput(1, 0) : indexOutput(1, -1);
	return {upper, lower};
}

EMatrix HydrogenFromFiles::cvv(double T, const EVector& speciesNv) const
{
	// Calculate the temperature in erg and in electron volt
	double kT = Constant::BOLTZMAN * T;
	double T_eV = kT * Constant::ERG_EV;

	EMatrix the_cvv = EMatrix::Zero(_numL, _numL);
	// Electron contributions (n-changing)
	double ne = speciesNv(_ine);
	if (ne > 0)
	{
		for (size_t i = 0; i < _numL; i++)
		{
			const HydrogenLevel& ini = _levelOrdering[i];
			for (size_t f = 0; f < _numL; f++)
			{
				const HydrogenLevel& fin = _levelOrdering[f];
				/* For downward transitions, calculate the collision rate, and
				   derive the rate for the corresponding upward transition too. */
				if (ini.e() > fin.e())
				{
					double UpsilonDown = eCollisionStrength(ini, fin, T_eV);
					double Cif = UpsilonDown * 8.6291e-6 / ini.g() /
					             sqrt(T) * ne;
					double Cfi = Cif * ini.g() / fin.g() *
					             exp((fin.e() - ini.e()) / kT);
					the_cvv(i, f) += Cif;
					the_cvv(f, i) += Cfi;
				}
				/* for upward transitions do nothing because we already covered them
				   above */
			}
		}
	}

	// For the l-resolved levels, get l-changing collision rates (through proton collisions)
	double np = speciesNv(_inp);
	if (np > 0)
	{
		for (int n = 0; n <= _resolvedUpTo; n++)
		{
			/* The formulae of PS64 are implemented separately. The result is indexed on
			   li, lf. */
			EMatrix qvv = PS64CollisionRateCoeff(n, T, np);
			// Fill in the collision rates for all combinations of li lf
			for (int li = 0; li < n; li++)
			{
				size_t i = indexOutput(n, li);
				for (int lf = 0; lf < n; lf++)
				{
					size_t f = indexOutput(n, lf);
					/* None of the previous contributions should have been
					   l-changing */
					assert(the_cvv(i, f) == 0);
					the_cvv(i, f) += qvv(li, lf) * np;
				}
			}
		}
	}
#ifdef SANITY
	if (!(the_cvv.array() >= 0).all())
		DEBUG("NEGATIVE COLLISION RATE");
#endif
	// Make all negative entries 0... FIXME investigate if a more better solution exists
	return the_cvv.array().max(0);
}

EMatrix HydrogenFromFiles::PS64CollisionRateCoeff(int n, double T, double np) const
{
	EMatrix q_li_lf_goingUp = EMatrix::Zero(n, n);
	EMatrix q_li_lf_goingDown = EMatrix::Zero(n, n);

	/* We will apply PS64 eq 43. Keeping eq 38 in mind, we can find the partial rates one by
	   one, applying detailed balance at each step. Note however that the results are different
	   depending on which side (l = 0 or l = n - 1) you start. Maybe we should calculate the
	   coefficients in both ways and then take the average. */

	// mu is the reduced mass of the system of the colliding particles
	constexpr double muOverm = Constant::PROTONMASS * Constant::HMASS_CGS /
	                           (Constant::PROTONMASS + Constant::HMASS_CGS) /
	                           Constant::ELECTRONMASS;

	const double qnlFactor = 9.93e-6 * sqrt(muOverm / T);
	const int n2 = n * n;

	auto ps64eq43 = [&](int l) -> double {
		// eq 44: Z is charge of the colliding particle, z that of the nucleus
		double D_nl = 6 * n2 * (n2 - l * l - l - 1);

		/* eq 45,46: take the smallest of the two, since Rc represents a cutoff value that
		   prevented divergence in the calculations of PS64 */
		size_t index = indexOutput(n, l);
		double tau2 = 1. / _totalAv(index) / _totalAv(index);
		double twoLog10Rc =
		                min(10.95 + log10(T * tau2 / muOverm), 1.68 + log10(T / np));

		// eq 43
		double q_nl = qnlFactor * D_nl *
		              (11.54 + log10(T / D_nl / muOverm) + twoLog10Rc);
		return max(0., q_nl);
	};

	double qDown = 0;
	// The value for q(l = n - 1, l = n - 2) is filled in in the last iteration
	for (int l = 0; l < n - 1; l++)
	{
		double qBoth = ps64eq43(l);
		double qUp = qBoth - qDown;
		q_li_lf_goingUp(l, l + 1) = qUp;

		// will also be used in loop for l+1:
		qDown = qUp * (2. * l + 1.) / (2. * l + 3.);
		q_li_lf_goingUp(l + 1, l) = qDown;
	}
	double qUp = 0;
	for (int l = n - 1; l > 0; l--)
	{
		double qBoth = ps64eq43(l);
		double qDown = qBoth - qUp;
		q_li_lf_goingDown(l, l - 1) = qDown;

		// will also be used in loop for l-1:
		qUp = qDown * (2. * l + 1.) / (2. * l - 1.);
		q_li_lf_goingDown(l - 1, l) = qUp;
	}
	const EMatrix result = (q_li_lf_goingUp + q_li_lf_goingDown) / 2.;
#ifdef SANITY
	assert((result.array() >= 0).all());
#endif
	return result.array().max(0);
}

EVector HydrogenFromFiles::sourcev(double T, const EVector& speciesNv) const
{
	/* TODO: find better recombination coefficients. */

	/* for now use the hardcoded implementation, but this needs to change (is copy paste from
	   HydrogenHardcoded). */
	double T4 = T / 1.e4;
	double alphaGround = 1.58e-13 * pow(T4, -0.53 - 0.17 * log(T4));
	double alpha2p = 5.36e-14 * pow(T4, -0.681 - 0.061 * log(T4));
	double alpha2s = 2.34e-14 * pow(T4, -0.537 - 0.019 * log(T4));

	// 2015-Raga (A13)
	double t = log10(T4);
	vector<double> logAlpha3poly = {-13.3377, -0.7161, -0.1435, -0.0386, 0.0077};
	vector<double> logAlpha4poly = {-13.5225, -0.7928, -0.1749, -0.0412, 0.0154};
	vector<double> logAlpha5poly = {-13.6820, -0.8629, -0.1957, -0.0375, 0.0199};

	double alpha3 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha3poly));
	double alpha4 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha4poly));
	double alpha5 = pow(10., TemplatedUtils::evaluatePolynomial(t, logAlpha5poly));

	// this hack should work for n up to 5
	Array unresolvedAlphav({alphaGround, alpha2p + alpha2s, alpha3, alpha4, alpha5});

	EVector alphav = EVector::Zero(_numL);

	// Now loop over all levels, and add the correct recombination coefficient
	for (size_t f = 0; f < _numL; f++)
	{
		const HydrogenLevel& final = _levelOrdering[f];
		int n = final.n();

		// Get one of the alphas calculated above
		double unresolvedAlpha = unresolvedAlphav[n - 1];

		// If it's the one for n = 2, use the resolved ones instead
		if (n == 2)
		{
			int l = final.l();
			if (l == 0)
				alphav[f] = alpha2s;
			else if (l == 1)
				alphav[f] = alpha2p;
			else if (l == -1)
				alphav[f] = unresolvedAlpha;
			else
				Error::runtime("invalid l value");
		}
		else
		{
			// weigh by multiplicity
			alphav[f] = unresolvedAlpha / (2 * n * n) * final.g();
		}
	}
	double ne = speciesNv(SpeciesIndex::ine());
	double np = speciesNv(SpeciesIndex::inp());
	return alphav * ne * np;
}

EVector HydrogenFromFiles::sinkv(double T, double n, const EVector& speciesNv) const
{
	/* The ionization rate calculation makes no distinction between the levels.  When
	   the upper level population is small, and its decay rate is large, the second term
	   doesn't really matter. Therefore, we choose the sink to be the same for each
	   level.  Moreover, total source = total sink so we want sink*n0 + sink*n1 = source
	   => sink = totalsource / n because n0/n + n1/n = 1. */
	double ne = speciesNv(_ine);
	double np = speciesNv(_inp);
	double totalSource = ne * np * Ionization::recombinationRateCoeff(T);
	double sink = totalSource / n; // Sink rate per (atom per cm3)
	return EVector::Constant(_numL, sink);
}

double HydrogenFromFiles::energy(int n, int l) const
{
	// Take an average over the j states
	double esum = 0;
	for (int twoJplus1 : twoJplus1range(l))
		esum += _chiantiLevelv[indexCHIANTI(n, l, twoJplus1)].e() * twoJplus1;
	return esum / (4 * l + 2);
}

double HydrogenFromFiles::energy(int n) const
{
	// Average over the l states
	double esum = 0;
	for (int l = 0; l < n; l++)
		esum += energy(n, l) * (2 * l + 1);
	return esum / (n * n);
}

double HydrogenFromFiles::einsteinA(int ni, int li, int nf, int lf) const
{
	if (ni < nf)
		return 0.;
	else
	{
		double Asum = 0;
		// sum over the final j states
		for (int twoJplus1f : twoJplus1range(lf))
		{
			int lIndex = indexCHIANTI(nf, lf, twoJplus1f);

			// average over the initial j states
			for (int twoJplus1i : twoJplus1range(li))
			{
				int uIndex = indexCHIANTI(ni, li, twoJplus1i);
				Asum += _chiantiAvv(uIndex, lIndex) * twoJplus1i;
			}
		}
		// divide by multiplicity of li state to get the average
		return Asum / (4 * li + 2);
	}
}

double HydrogenFromFiles::einsteinA(int ni, int nf, int lf) const
{
	// average over the initial l states
	double Asum = 0;
	for (int li = 0; li < ni; li++)
		Asum += einsteinA(ni, li, nf, lf) * (2 * li + 1);
	return Asum / (ni * ni);
}

double HydrogenFromFiles::einsteinA(int ni, int nf) const
{
	// sum over the final l states
	double Asum = 0;
	for (int lf = 0; lf < nf; lf++)
		Asum += einsteinA(ni, nf, lf);
	return Asum;
}

double HydrogenFromFiles::einsteinA(const HydrogenLevel& initial,
                                    const HydrogenLevel& final) const
{
	// No output for upward transitions
	if (initial.e() < final.e())
		return 0;
	else
	{
		// Resolved-resolved
		if (initial.nlResolved() && final.nlResolved())
			return einsteinA(initial.n(), initial.l(), final.n(), final.l());
		// Collapsed-resolved
		else if (initial.nCollapsed() && final.nlResolved())
			return einsteinA(initial.n(), final.n(), final.l());
		/* Resolved-collapsed (should not be called, and if it is, the collapsed-collapsed
		   result should be zero). */
		else if (initial.nlResolved() && final.nCollapsed())
		{
			double a = einsteinA(initial.n(), final.n());
			assert(a == 0);
			return a;
		}
		// Collapsed-collapsed
		else
			return einsteinA(initial.n(), final.n());
	}
}

double HydrogenFromFiles::eCollisionStrength(int ni, int li, int nf, int lf, double T_eV) const
{
	/* Find the requested transition in the map that translates n,l to the index in the Anderson
	   file. */
	auto uAndersonIndexIt = _nlToAndersonIndexm.find({ni, li});
	auto lAndersonIndexIt = _nlToAndersonIndexm.find({nf, lf});
	auto indexMapEnd = _nlToAndersonIndexm.end();
	/* When either of the levels is not in the map, that means there are no transitions
	   involving this level in the data set. The result is 0. */
	if (uAndersonIndexIt == indexMapEnd || lAndersonIndexIt == indexMapEnd)
		return 0;

	int uAndersonIndex = uAndersonIndexIt->second;
	int lAndersonIndex = lAndersonIndexIt->second;
	if (uAndersonIndex <= lAndersonIndex)
		Error::runtime("This function should only be used for downward collisional "
		               "transitions");

	// Try to find this transition in the data set.
	auto UpsilonvIt = _andersonUpsilonvm.find({uAndersonIndex, lAndersonIndex});
	// If there is no entry, the result is 0.
	if (UpsilonvIt == _andersonUpsilonvm.end())
		return 0;

	/* If data is available for this transition, interpolate the Upsilon(T)-array (at the given
	   temperature (in eV). */
	double Upsilon_ip = TemplatedUtils::evaluateLinInterpf<double, Array, Array>(
	                T_eV, _andersonTempv, UpsilonvIt->second);
	return Upsilon_ip;
}

double HydrogenFromFiles::eCollisionStrength(int ni, int nf, int lf, double T_eV) const
{
	// Collision strengths must be summed over all initial states
	double Upsilonsum = 0;
	for (int li = 0; li < ni; li++)
		Upsilonsum += eCollisionStrength(ni, li, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenFromFiles::eCollisionStrength(int ni, int nf, double T_eV) const
{
	// Also sum over all final states
	double Upsilonsum = 0;
	for (int lf = 0; lf < nf; lf++)
		Upsilonsum += eCollisionStrength(ni, nf, lf, T_eV);
	return Upsilonsum;
}

double HydrogenFromFiles::eCollisionStrength(const HydrogenLevel& initial,
                                             const HydrogenLevel& final, double T_eV) const
{
	// No output for upward transitions
	if (initial.e() < final.e())
		return 0;
	else
	{
		// Resolved-resolved
		if (initial.nlResolved() && final.nlResolved())
			return eCollisionStrength(initial.n(), initial.l(), final.n(),
			                          final.l(), T_eV);
		// Collapsed-resolved
		else if (initial.nCollapsed() && final.nlResolved())
			return eCollisionStrength(initial.n(), final.n(), final.l(), T_eV);
		/* Resolved-collapsed (should not be called as this does not exist for downward
		   transitions, and if it is, the collapsed-collapsed result should be zero. */
		else if (initial.nlResolved() && final.nCollapsed())
		{
			double eCollStr = eCollisionStrength(initial.n(), final.n(), T_eV);
			assert(eCollStr == 0);
			return eCollStr;
		}
		// Collapsed-collapsed
		else
			return eCollisionStrength(initial.n(), final.n(), T_eV);
	}
}
