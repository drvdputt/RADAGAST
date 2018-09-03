#include "GasInterfaceImpl.h"
#include "ChemistrySolver.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "GasGrainInteraction.h"
#include "GasStruct.h"
#include "GrainType.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "SimpleHydrogenNetwork.h"
#include "SpecialFunctions.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"
#include "Testing.h"

using namespace std;

constexpr int MAXCHEMISTRYITERATIONS{100};
constexpr bool GASGRAINCOOL{true};

GasInterfaceImpl::GasInterfaceImpl(unique_ptr<HydrogenLevels> atomModel,
                                   unique_ptr<H2Levels> molecularModel)
                : _atomicLevels(move(atomModel)), _molecular(move(molecularModel)),
                  _freeBound(make_unique<FreeBound>()), _freeFree(make_unique<FreeFree>())
{
	if (_molecular)
		_chemSolver = make_unique<ChemistrySolver>(
		                make_unique<SimpleHydrogenNetwork>());

	_ine = SpeciesIndex::index("e-");
	_inp = SpeciesIndex::index("H+");
	_inH = SpeciesIndex::index("H");
	_inH2 = SpeciesIndex::index("H2");
}

GasInterfaceImpl::~GasInterfaceImpl() = default;

void GasInterfaceImpl::solveInitialGuess(GasModule::GasState& gs, double n, double T,
                                         const GasModule::GrainInterface& gi,
                                         const Array& oFrequencyv,
                                         const Array& eFrequencyv) const
{
	Array isrfGuess(oFrequencyv.size());
	for (size_t iFreq = 0; iFreq < oFrequencyv.size(); iFreq++)
	{
		double freq = oFrequencyv[iFreq];
		isrfGuess[iFreq] = SpecialFunctions::planck(freq, T);
	}
	solveBalance(gs, n, T, Spectrum(oFrequencyv, isrfGuess), gi, oFrequencyv, eFrequencyv);
}

void GasInterfaceImpl::solveBalance(GasModule::GasState& gs, double n, double Tinit,
                                    const Spectrum& specificIntensity,
                                    const GasModule::GrainInterface& gi,
                                    const Array& oFrequencyv, const Array& eFrequencyv) const
{
#ifndef SILENT
	double isrf = TemplatedUtils::integrate<double>(specificIntensity.frequencyv(),
	                                                specificIntensity.valuev());

	DEBUG("Solving balance under isrf of " << isrf << " erg / s / cm2 / sr = ");
	DEBUG(isrf / Constant::LIGHT * Constant::FPI / Constant::HABING << " Habing" << endl);
#endif
	Solution s;

	if (n <= 0)
		s = calculateDensities(0, 0, specificIntensity, gi);

	else
	{
		const double Tmax = 1000000.;
		const double Tmin = 10;
		const double logTmax = log10(Tmax);
		const double logTmin = log10(Tmin);
		const double logTtolerance = 1.e-3;

		double logTinit = log10(Tinit);

		/* Lambda function that will be used by the search algorithm. The state of the
		   system will be updated every time the algorithm calls this function. The
		   return value indicates whether the temperature should increase (there is net
		   heating so we need a higher temperature leading to more cooling) or decrease
		   (there is net cooling so we need a lower temperature leading to less
		   cooling). */
		int counter = 0;
		const Solution* previous = nullptr;
		function<int(double)> evaluateThermalBalance = [&](double logT) -> int {
			counter++;
			s = calculateDensities(n, pow(10., logT), specificIntensity, gi,
			                       previous);
			previous = &s;
			double netPowerIn = heating(s, gi) - cooling(s);
			DEBUG("Cycle " << counter << ": logT = " << logT
			               << "; netHeating = " << netPowerIn << endl
			               << endl);
			return (netPowerIn > 0) - (netPowerIn < 0);
		};

		double logTfinal = TemplatedUtils::binaryIntervalSearch<double>(
		                evaluateThermalBalance, logTinit, logTtolerance, logTmax,
		                logTmin);

		// Evaluate the densities for one last time, using the final temperature.
		s = calculateDensities(n, pow(10., logTfinal), specificIntensity, gi, previous);
	}

	const Array& emv = emissivityv(s, eFrequencyv);
	const Array& opv = opacityv(s, oFrequencyv);
	const Array& scv = scatteringOpacityv(s, oFrequencyv);

#ifdef SANITY
	for (size_t iFreq = 0; iFreq < _frequencyv.size(); iFreq++)
	{
		if (emv[iFreq] < 0 || opv[iFreq] < 0 || scv[iFreq] < 0)
		{
			Error::runtime("GasModule: negative value in one of the optical "
			               "properties");
		}
	}
#endif
	// Put the relevant data into the gas state
	Array densityv(s.speciesNv.data(), s.speciesNv.size());
	// Derive this again, just for diagnostics
	double h2form = GasGrain::surfaceH2FormationRateCoeff(gi, s.T) * s.speciesNv[_inH];
	double grainHeat = grainHeating(s, gi);
	gs = GasModule::GasState(emv, opv, scv, s.T, densityv, h2form, grainHeat);
}

GasInterfaceImpl::Solution
GasInterfaceImpl::calculateDensities(double nHtotal, double T,
                                     const Spectrum& specificIntensity,
                                     const GasModule::GrainInterface& gi,
                                     const GasInterfaceImpl::Solution* previous) const
{
	Solution s;

	// We can already fill these in.
	s.T = T;
	s.specificIntensity = specificIntensity;

	if (nHtotal <= 0)
	{
		s.speciesNv = EVector::Zero(SpeciesIndex::size());
		s.HSolution = _atomicLevels->solveZero(T);
		if (_molecular)
			s.H2Solution = _molecular->solveZero(T);
		return s;
	}

	DEBUG("Calculating densities for T = " << T << "K" << endl);

	// Decide how to do the initial guess; manual, or using the previous state.
	bool manualGuess = true;

	// If a previous solution is available (pointer is nonzero)
	if (previous)
	{
		double fT = T / previous->T;
		if (fT > 0.2 && fT < 5)
		{
			DEBUG("Using previous speciesNv as initial guess" << std::endl);
			s.speciesNv = previous->speciesNv;
			manualGuess = false;
		}
		// If the difference in temperature is larger than a certain factor, do a manual
		// guess anyway.
	}
	if (manualGuess)
	{
		double iniNH2 = nHtotal / 4;
		double iniAtomAndIon = nHtotal / 2;
		double guessF = Ionization::solveBalance(iniAtomAndIon, T, specificIntensity);
		s.speciesNv = EVector(SpeciesIndex::size());
		s.speciesNv(_ine) = guessF * iniAtomAndIon;
		s.speciesNv(_inp) = guessF * iniAtomAndIon;
		s.speciesNv(_inH) = (1 - guessF) * iniAtomAndIon;
		s.speciesNv(_inH2) = iniNH2;
	}

	/* Lambda function, because it is only needed in this scope. The [&] passes the current
	   scope by reference, so the lambda function can modify s. */
	// Package some gas parameters
	GasStruct gas(T, s.speciesNv);

	// Initial guess for the H2 solution
	if (_molecular)
	{
		if (manualGuess)
		{
			DEBUG("Using LTE as initial guess for H2" << endl);
			gas._h2Levelv = s.speciesNv(_inH2) *
			                _molecular->solveBoltzmanEquations(gas._T);
		}
		else
			gas._h2Levelv = previous->H2Solution.nv;
	}

	auto solveLevelBalances = [&]() {
		double nH = s.speciesNv(_inH);
		EVector Hsourcev = _atomicLevels->sourcev(gas);
		EVector Hsinkv = _atomicLevels->sinkv(gas);

		DEBUG("Solving levels nH = " << nH << endl);
		s.HSolution = _atomicLevels->solveBalance(nH, specificIntensity, Hsourcev,
		                                          Hsinkv, gas);
		if (_molecular)
		{
			double nH2 = s.speciesNv(_inH2);
			// TODO: the source term should contain the 'formation pumping'
			// contributions. When H2 is formed on grain surfaces, it can be
			// released from the grain in an excited state. The simplest way is
			// assuming a fixed distribution. In that case, the source vector is
			// this distribution scaled with the total grain H2 formation rate.
			EVector H2sourcev = EVector::Zero(_molecular->numLv());
			EVector H2sinkv = _molecular->dissociationSinkv(specificIntensity);
			DEBUG("Solving levels nH2 = " << nH2 << endl);
			s.H2Solution = _molecular->solveBalance(nH2, specificIntensity,
			                                        H2sourcev, H2sinkv, gas);

			// Save this, to be used as an initial condition later
			gas._h2Levelv = s.H2Solution.nv;
		}
	};

	/* The main loop: e
         v---------------------------------------------------<
	 > LEVELS SOLUTION, CHEMISTRY SOLUTION -> CHEM RATES |
                                                             |
	   CHEM RATES -> CHEMISTRY SOLUTION                  |
                                                             |
	   CHEMISTRY SOLUTION -> LEVEL SOURCE AND SINK RATES |
                                                             |
	   LEVEL SOURCE AND SINK RATES -> LEVEL SOLUTION ----^
	*/

	int counter = 0;
	bool stopCriterion = false;
	while (!stopCriterion)
	{
		/* Use the current guess for the chemical abundances to calculate our first set
		   of level populations. */

		// CHEMISTRY SOLUTION -> SOURCE AND SINK RATES -> LEVEL SOLUTION
		solveLevelBalances();

		EVector previousAbundancev = s.speciesNv;

		if (previousAbundancev.array().isNaN().any())
			Error::runtime("Nan in chemistry solution!");

		// When including H2
		if (_molecular)
		{
			// LEVELS AND CHEMISTRY SOLUTIONS -> CHEM RATES
			double kFormH2 = GasGrain::surfaceH2FormationRateCoeff(gi, T);
			double kDissH2Levels = _molecular->dissociationRate(
			                s.H2Solution, s.specificIntensity);

			DEBUG("Formation rate per H " << kFormH2 << endl);
			DEBUG("Dissociation rate per H2 " << kDissH2Levels << endl);
			if (kDissH2Levels < 0)
				Error::runtime("negative dissociation rate!");

			// CHEM RATES -> CHEMISTRY SOLUTION
			EVector reactionRates = _chemSolver->chemicalNetwork()->rateCoeffv(
			                T, specificIntensity, kDissH2Levels, kFormH2);
			s.speciesNv = _chemSolver->solveBalance(reactionRates, s.speciesNv);

			/* TODO: Add effect of grain charging to chemical network. I think it
			   might be possible to do this by imposing a conservation equation for
			   the number of electrons: ne + nH + nH2 = (ne + nH + nH2)_0 + <Cg>*ng.
			   The average grain charge <Gg> should be updated together with the
			   rates I guess? Another option would be to include the grain charge
			   rates into the network as extra reactions. The production vector
			   would be (1 0 0 0) while the reactant vector would be zero (the
			   grains don't disappear when they lose an electron) Grain
			   recombination / charge exchange reaction could also be added. I need
			   to think about whether the 'disappearing' particles will cause
			   problems when coupled with conservation equations. */
		}
		// When ignoring H2
		else
		{
			// DIRECT SOLUTION (IONIZATION BALANCE ONLY, NO MOLECULES)

			// Just solve the ionization balance in the nebular approximation.
			double f = Ionization::solveBalance(nHtotal, T, specificIntensity);
			DEBUG("Ionized fraction = " << f << endl);

			// Neutral fraction
			s.speciesNv[_inH] = nHtotal * (1 - f);
			// Ionized fraction
			s.speciesNv[_inp] = nHtotal * f;
			// Electron density is simply equal to proton density
			s.speciesNv[_ine] = s.speciesNv[_inp];
			s.speciesNv[_inH2] = 0;
		}

		// CONVERGENCE CHECK. An abundance has converged if it changes by less than 1%
		EArray changev = s.speciesNv - previousAbundancev;
		// Or if it is negligible compared to the norm (or sum maybe?)
		double norm = s.speciesNv.norm();
		Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
		                changev.abs() <= 0.01 * previousAbundancev.array() ||
		                s.speciesNv.array() < 1.e-99 * norm;
		counter++;
		DEBUG("Chemistry: " << counter << endl
		                    << s.speciesNv << endl
		                    << "convergence: \n"
		                    << convergedv << endl);

		// Currently, the implementation without molecules does not need iteration.
		stopCriterion = !_molecular || convergedv.all() ||
		                counter > MAXCHEMISTRYITERATIONS;
	}
	return s;
}

Array GasInterfaceImpl::emissivityv(const Solution& s, const Array& eFrequencyv) const
{
	Array lineEmv(eFrequencyv.size());
	lineEmv = _atomicLevels->emissivityv(s.HSolution, eFrequencyv);
	if (_molecular)
		lineEmv += _molecular->emissivityv(s.H2Solution, eFrequencyv);

	Array contEmCoeffv(eFrequencyv.size());
	_freeBound->addEmissionCoefficientv(s.T, eFrequencyv, contEmCoeffv);
	_freeFree->addEmissionCoefficientv(s.T, eFrequencyv, contEmCoeffv);

	return lineEmv + (np_ne(s) / Constant::FPI) * contEmCoeffv;
}

Array GasInterfaceImpl::opacityv(const Solution& s, const Array& oFrequencyv) const
{
	size_t numFreq = oFrequencyv.size();
	Array lineOp(numFreq);
	lineOp += _atomicLevels->opacityv(s.HSolution, oFrequencyv);
	if (_molecular)
		lineOp += _molecular->opacityv(s.H2Solution, oFrequencyv);

	Array contOpCoeffv(numFreq);
	_freeFree->addOpacityCoefficientv(s.T, oFrequencyv, contOpCoeffv);

	double npne = np_ne(s);
	double nH0 = nAtomic(s);
	Array totalOp(numFreq);
	for (size_t iFreq = 0; iFreq < numFreq; iFreq++)
	{
		double ionizOp_iFreq = nH0 * Ionization::crossSection(oFrequencyv[iFreq]);
		totalOp[iFreq] = ionizOp_iFreq + npne * contOpCoeffv[iFreq] + lineOp[iFreq];
#ifdef SANITY
		if (totalOp[iFreq] < 0)
			cout << "Negative opacity!";
#endif
	}
	return totalOp;
}

Array GasInterfaceImpl::scatteringOpacityv(const Solution& s, const Array& oFrequencyv) const
{
	return Array(oFrequencyv.size());
}

double GasInterfaceImpl::cooling(const Solution& s) const
{
	double lineCool = lineCooling(s);
	double contCool = continuumCooling(s);
	DEBUG("cooling: line / cont = " << lineCool << " / " << contCool << endl);
	return lineCool + contCool;
}

double GasInterfaceImpl::heating(const Solution& s, const GasModule::GrainInterface& g) const
{
	return heating(s) + grainHeating(s, g);
}

double GasInterfaceImpl::heating(const Solution& s) const
{
	double lineHeat = lineHeating(s);
	double contHeat = continuumHeating(s);
	DEBUG("heating: line / cont = " << lineHeat << " / " << contHeat << endl);
	return lineHeat + contHeat;
}

double GasInterfaceImpl::grainHeating(const Solution& s,
                                      const GasModule::GrainInterface& g) const
{
	double grainPhotoelectricHeating{0};
	// Specify the environment parameters
	double ne = s.speciesNv[_ine];
	double np = s.speciesNv[_inp];
	GrainPhotoelectricEffect::Environment env(
	                s.specificIntensity, s.T, ne, np, {-1, 1}, {ne, np},
	                {Constant::ELECTRONMASS, Constant::PROTONMASS});
	size_t numPop = g.numPopulations();
	for (size_t i = 0; i < numPop; i++)
	{
		const GasModule::GrainInterface::Population* pop = g.population(i);
		const GrainType* type = pop->type();
		if (type->heatingAvailable())
		{
			/* Choose the correct parameters for the photoelectric heating
			   calculation based on the type (a.k.a. composition) of the
			   Population. */
			GrainPhotoelectricEffect gpe(*type);
			grainPhotoelectricHeating += gpe.heatingRate(env, *pop);
		}
	}

	// For gas-grain collisional cooling, calculate the average dust temperature, weighted
	// by cross section. This will be the expected value of the temperatures a particle
	// 'experiences' when colliding with a grain.
	double Tsum = 0;
	double weight = 0;
	for (size_t i = 0; i < numPop; i++)
	{
		const GasModule::GrainInterface::Population* pop = g.population(i);
		for (size_t m = 0; m < pop->numSizes(); m++)
		{
			double a = pop->size(m);
			double nd = pop->density(m);
			double w = Constant::PI * a * a * nd;
			Tsum += w * pop->temperature(m);
			weight += w;
		}
	}
	double Tdust = Tsum / weight;
	double Tgas = s.T;
	double nH = s.speciesNv[_inH];
	double nH2 = s.speciesNv[_inH2];
	double gasGrainCooling =
	                GASGRAINCOOL ? GasGrain::simpleGasGrainCool(Tdust, Tgas, nH, nH2) : 0;

	return grainPhotoelectricHeating - gasGrainCooling;
}

double GasInterfaceImpl::lineCooling(const Solution& s) const
{
	double result = _atomicLevels->cooling(s.HSolution);
	if (_molecular)
		result += _molecular->cooling(s.H2Solution);
	return result;
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	double result = _atomicLevels->heating(s.HSolution);
	if (_molecular)
		result += _molecular->heating(s.H2Solution);
	return result;
}

double GasInterfaceImpl::continuumCooling(const Solution& s) const
{
	double result = _freeFree->cooling(np_ne(s), s.T) +
	                Ionization::cooling(s.speciesNv(_inH), s.speciesNv(_inp),
	                                    s.speciesNv(_ine), s.T);
	if (_molecular)
		result += _molecular->dissociationCooling(s.H2Solution);
	return result;
}

double GasInterfaceImpl::continuumHeating(const Solution& s) const
{
	double result = _freeFree->heating(np_ne(s), s.T, s.specificIntensity);
	result += Ionization::heating(s.speciesNv(_inp), s.speciesNv(_ine), s.T,
	                              s.specificIntensity);
	if (_molecular)
	{
		double dissheat = _molecular->dissociationHeating(s.H2Solution);
		result += dissheat;
	}
	return result;
}
