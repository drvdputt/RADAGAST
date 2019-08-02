#include "GasInterfaceImpl.h"
#include "ChemistrySolver.h"
#include "Constants.h"
#include "DebugMacros.h"
#include "FreeBound.h"
#include "FreeFree.h"
#include "GasDiagnostics.h"
#include "GasGrainInteraction.h"
#include "GasStruct.h"
#include "GrainType.h"
#include "H2FromFiles.h"
#include "H2Levels.h"
#include "HydrogenFromFiles.h"
#include "HydrogenLevels.h"
#include "IOTools.h"
#include "IonizationBalance.h"
#include "LevelSolver.h"
#include "SimpleHydrogenNetwork.h"
#include "SpecialFunctions.h"
#include "SpeciesIndex.h"
#include "TemplatedUtils.h"
#include "Testing.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

using namespace std;

constexpr int MAXCHEMISTRYITERATIONS{25};
constexpr bool GASGRAINCOOL{false};

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

GasInterfaceImpl::Solution
GasInterfaceImpl::solveInitialGuess(double n, double T,
                                    const GasModule::GrainInterface& gi) const
{
	Array iFrequencyv = Testing::defaultFrequencyv(1000);
	Array blackbodyv(iFrequencyv.size());
	for (size_t i = 0; i < iFrequencyv.size(); i++)
		blackbodyv[i] = SpecialFunctions::planck(iFrequencyv[i], T);

	Solution s = solveTemperature(n, T, Spectrum(iFrequencyv, blackbodyv), gi);
	return s;
}

struct heating_f_params
{
	const GasInterfaceImpl* gasInterfacePimpl;
	double n;
	const Spectrum* specificIntensity;
	const GasModule::GrainInterface* grainInterface;
	GasInterfaceImpl::Solution* solution_storage;
	bool use_previous_solution;
};

/* Function that will be used by the GSL search algorithm to find the equilibrium temperature.
   The state of the system will be updated every time the algorithm calls this function (and
   stored via the pointer suppied in the heating_f_params struct). */
double heating_f(double logT, void* params)
{
	auto* p = static_cast<struct heating_f_params*>(params);

	// Refresh the solution stored somewhere, with optimization based on the current
	// solution
	GasInterfaceImpl::Solution& s = *p->solution_storage;
	const GasInterfaceImpl::Solution* previous =
	                p->use_previous_solution ? p->solution_storage : nullptr;

	s = p->gasInterfacePimpl->solveDensities(p->n, pow(10., logT), *p->specificIntensity,
	                                         *p->grainInterface, previous);

	double heat = p->gasInterfacePimpl->heating(s, *p->grainInterface);
	double cool = p->gasInterfacePimpl->cooling(s);
	return heat - cool;
}

GasInterfaceImpl::Solution
GasInterfaceImpl::solveTemperature(double n, double /*unused Tinit*/,
                                   const Spectrum& specificIntensity,
                                   const GasModule::GrainInterface& gri) const
{
	Solution s;

	if (n <= 0)
		s = solveDensities(0, 0, specificIntensity, gri);

	else
	{
		const double Tmax = 1000000.;
		const double Tmin = 10;
		double logTmax = log10(Tmax);
		const double logTmin = log10(Tmin);
		const double logTtolerance = 1.e-3;

		const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
		gsl_root_fsolver* solver = gsl_root_fsolver_alloc(T);

		gsl_function F;
		struct heating_f_params p = {this, n, &specificIntensity, &gri, &s, false};
		F.function = &heating_f;
		F.params = &p;

		// For very high temperatures, the heating becomes positive again, and the GSL
		// function will error out. It needs to to have different signs for the points
		// bracketing the root.
		double reduceMax = 1.;
		while (heating_f(logTmax, &p) > 0)
		{
			while (logTmax - reduceMax < logTmin)
				reduceMax /= 2;

			logTmax -= reduceMax;
		}

		// Iterate once to initialize the solution
		gsl_root_fsolver_set(solver, &F, logTmin, logTmax);
		gsl_root_fsolver_iterate(solver);

		// Then, start using the current solution as an initial guess of the next one
		p.use_previous_solution = true;

		gsl_root_fsolver_set(solver, &F, logTmin, logTmax);
		int test_interval = GSL_CONTINUE;
		int counter = 0;
		while (test_interval != GSL_SUCCESS)
		{
			gsl_root_fsolver_iterate(solver);
			double lower = gsl_root_fsolver_x_lower(solver);
			double upper = gsl_root_fsolver_x_upper(solver);
			test_interval = gsl_root_test_interval(lower, upper, logTtolerance,
			                                       logTtolerance);
			counter++;
		}
		gsl_root_fsolver_free(solver);
		// Remember that s gets updated through the pointer given to the function
		// parameters
	}
	return s;
}

GasModule::GasState GasInterfaceImpl::makeGasStateFromSolution(const Solution& s,
                                                               const Array& oFrequencyv,
                                                               const Array& eFrequencyv) const

{
	Array emv, opv, scv;
	if (eFrequencyv.size() > 2)
		emv = emissivityv(s, eFrequencyv);
	if (oFrequencyv.size() > 2)
	{
		opv = opacityv(s, oFrequencyv);
		scv = scatteringOpacityv(s, oFrequencyv);
	}
	Array densityv(s.speciesNv.data(), s.speciesNv.size());
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
	return {emv, opv, scv, s.T, densityv};
}

void GasInterfaceImpl::fillGasDiagnosticsFromSolution(const Solution& s,
                                                      const GasModule::GrainInterface& gri,
                                                      GasDiagnostics* gd) const
{
	if (!gd)
		Error::runtime("GasDiagnostics is nullptr!");

	gd->setPhotoelectricHeating(Array({grainHeating(s, gri)}));

	double h2form = GasGrain::surfaceH2FormationRateCoeff(gri, s.T);
	double h2dissoc = _molecular ? _molecular->dissociationRate(s.H2Solution,
	                                                            s.specificIntensity)
	                             : 0;
	double hphotoion = Ionization::photoRateCoeff(s.specificIntensity);
	double hcolion = Ionization::collisionalRateCoeff(s.T);
	double hrec = Ionization::recombinationRateCoeff(s.T);

	gd->setReactionNames({"h2form", "h2dissoc", "hphotoion", "hcolion", "hrec"});
	gd->setReactionRates({h2form, h2dissoc, hphotoion, hcolion, hrec});

	if (gd->saveLevelPopulations())
	{
		gd->setHPopulations(Array(s.HSolution.nv.data(), s.HSolution.nv.size()));
		gd->setH2Populations(Array(s.H2Solution.nv.data(), s.H2Solution.nv.size()));
	}

	gd->setHeating("H ion", Ionization::heating(s.speciesNv(_inp), s.speciesNv(_ine), s.T,
	                                            s.specificIntensity));
	gd->setHeating("H deexc", _atomicLevels->heating(s.HSolution));
	gd->setHeating("H2 deexc", _molecular ? _molecular->heating(s.H2Solution) : 0);
	gd->setHeating("H2 dissoc",
	               _molecular ? _molecular->dissociationHeating(s.H2Solution,
	                                                            s.specificIntensity)
	                          : 0);
}

GasInterfaceImpl::Solution
GasInterfaceImpl::solveDensities(double nHtotal, double T, const Spectrum& specificIntensity,
                                 const GasModule::GrainInterface& gi,
                                 const GasInterfaceImpl::Solution* previous,
                                 double h2FormationOverride) const
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

	// If a previous solution is available (pointer is nonzero), then check if we want to
	// use it
	if (previous)
	{
		double fT = T / previous->T;
		if (fT > 0.75 && fT < 1.5)
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
		double ionFrac = Ionization::solveBalance(nHtotal, T, specificIntensity);
		double molFrac = 0.1;
		s.speciesNv = EVector(SpeciesIndex::size());
		s.speciesNv(_ine) = ionFrac * nHtotal;
		s.speciesNv(_inp) = s.speciesNv(_ine);
		s.speciesNv(_inH) = (1 - molFrac) * (1 - ionFrac) * nHtotal;
		s.speciesNv(_inH2) = molFrac * (1 - ionFrac) * nHtotal / 2;
	}
	else
	{
		// make sure that nHtotal and the initial guess are consistent
		double nNeutralOfSpeciesNv = s.speciesNv(_inH) + 2 * s.speciesNv(_inH2);
		double nHtotalOfSpeciesNv = s.speciesNv(_inp) + nNeutralOfSpeciesNv;
		double ionFrac = s.speciesNv(_inp) / nHtotalOfSpeciesNv;
		double molFrac = 2 * s.speciesNv(_inH2) / nNeutralOfSpeciesNv;
		s.speciesNv(_ine) = ionFrac * nHtotal;
		s.speciesNv(_inp) = s.speciesNv(_ine);
		s.speciesNv(_inH) = (1 - molFrac) * (1 - ionFrac) * nHtotal;
		s.speciesNv(_inH2) = molFrac * (1 - ionFrac) * nHtotal / 2;
	}

	// Package some gas parameters
	GasStruct gas(T, s.speciesNv);

	// Formation rate on grain surfaces. We declare it here so that it can be passed to the
	// level population calculation too (needed for H2 formation pumping).
	double kFormH2 = 0;

	// Initial guess for the H2 solution (if no initial guess is provided, the solver will
	// make its own)
	if (_molecular && !manualGuess)
		gas._h2Levelv = nHtotal * previous->H2Solution.nv /
		                previous->H2Solution.nv.sum();

	auto solveLevelBalances = [&]() {
		double nH = s.speciesNv(_inH);
		if (nH <= 0)
			s.HSolution = _atomicLevels->solveZero(T);
		else
		{
			s.HSolution.n = nH;
			s.HSolution.T = T;
			EMatrix Htransitionvv = _atomicLevels->totalTransitionRatesvv(
			                specificIntensity, gas, &s.HSolution.cvv);
			EVector Hsourcev = _atomicLevels->sourcev(gas);
			EVector Hsinkv = _atomicLevels->sinkv(gas);
			DEBUG("Solving levels nH = " << nH << endl);
			s.HSolution.nv = LevelSolver::statisticalEquilibrium(nH, Htransitionvv,
			                                                     Hsourcev, Hsinkv);
			// s.HSolution = _atomicLevels->solveLTE(nH, gas);
		}

		if (_molecular)
		{
			double nH2 = s.speciesNv(_inH2);
			DEBUG("Solving levels nH2 = " << nH2 << endl);
			// gas is passed here to make an initial guess of the levels based on
			// the gas properties
			s.H2Solution = _molecular->customSolution(nH2, gas, specificIntensity,
			                                          nH * kFormH2);

			// the solution is also saved in the gas struct, so it can be used as an
			// initial guess instead of just guessing based on the temperature
			gas._h2Levelv = s.H2Solution.nv;
			// s.H2Solution = _molecular->solveLTE(nH2, gas);
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
	EVector previousAbundancev = s.speciesNv;
	double previousHeating = 0;
	double previousCooling = 0;
	while (!stopCriterion)
	{
		// CHEMISTRY SOLUTION -> SOURCE AND SINK RATES -> LEVEL SOLUTION
		solveLevelBalances();

		if (previousAbundancev.array().isNaN().any())
			Error::runtime("Nan in chemistry solution!");

		// When including H2
		if (_molecular)
		{
			// LEVELS AND CHEMISTRY SOLUTIONS -> CHEM RATES
			kFormH2 = GasGrain::surfaceH2FormationRateCoeff(gi, T);
			if (h2FormationOverride >= 0)
				kFormH2 = h2FormationOverride;

			double kDissH2Levels = _molecular->dissociationRate(
			                s.H2Solution, s.specificIntensity);

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
		                changev.abs() <= 1e-3 * previousAbundancev.array() ||
		                s.speciesNv.array() < 1.e-99 * norm;

		double newHeating = heating(s);
		double newCooling = cooling(s);
		bool heatingConverged = TemplatedUtils::equalWithinTolerance(
		                newHeating, previousHeating, 1e-2);
		bool coolingConverged = TemplatedUtils::equalWithinTolerance(
		                newCooling, previousCooling, 1e-2);
		counter++;
		DEBUG("Chemistry: " << counter << '\n'
		                    << s.speciesNv << '\n'
		                    << "convergence: \n"
		                    << convergedv << '\n');
		DEBUG("New heat: " << newHeating << " previous: " << previousHeating << '\n');
		DEBUG("New cool: " << newCooling << " previous: " << previousCooling << '\n');

		previousAbundancev = s.speciesNv;
		previousHeating = newHeating;
		previousCooling = newCooling;

		bool allQuantitiesConverged =
		                convergedv.all() && heatingConverged && coolingConverged;

		// Currently, the implementation without molecules does not need iteration.
		stopCriterion = !_molecular || allQuantitiesConverged ||
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

Array GasInterfaceImpl::scatteringOpacityv(const Solution&, const Array& oFrequencyv) const
{
	return Array(oFrequencyv.size());
}

Array GasInterfaceImpl::radiativeRecombinationEmissivityv(const Solution& s,
                                                          const Array& eFrequencyv) const
{
	Array rrEmv(eFrequencyv.size());
	_freeBound->addEmissionCoefficientv(s.T, eFrequencyv, rrEmv);
	return rrEmv;
}

Array GasInterfaceImpl::freeFreeEmissivityv(const Solution& s, const Array& eFrequencyv) const
{
	Array ffEmv(eFrequencyv.size());
	_freeFree->addEmissionCoefficientv(s.T, eFrequencyv, ffEmv);
	return ffEmv;
}

Array GasInterfaceImpl::lineEmissivityv(const Solution& s, const Array& eFrequencyv) const
{
	Array lineEmv(eFrequencyv.size());
	lineEmv = _atomicLevels->emissivityv(s.HSolution, eFrequencyv);
	if (_molecular)
		lineEmv += _molecular->emissivityv(s.H2Solution, eFrequencyv);
	return lineEmv;
}

Array GasInterfaceImpl::ionizationOpacityv(const Solution& s, const Array& oFrequencyv) const
{
	Array cs(oFrequencyv.size());
	for (size_t iFreq = 0; iFreq < oFrequencyv.size(); iFreq++)
		cs[iFreq] = s.speciesNv(SpeciesIndex::inH()) *
		            Ionization::crossSection(oFrequencyv[iFreq]);
	return cs;
}

// Array& GasInterfaceImpl::dissociationOpacityv(const Solution& s, const Array& oFrequencyv) const {
// 	Array cs(eFrequencyv.size());

// }

double GasInterfaceImpl::cooling(const Solution& s) const
{
	double lineCool = lineCooling(s);
	double contCool = continuumCooling(s);
	return lineCool + contCool;
}

double GasInterfaceImpl::heating(const Solution& s, const GasModule::GrainInterface& g) const
{
	double gasHeat = heating(s);
	double grainHeat = grainHeating(s, g);
	return gasHeat + grainHeat;
}

double GasInterfaceImpl::heating(const Solution& s) const
{
	double lineHeat = lineHeating(s);
	double contHeat = continuumHeating(s);
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

	double gasGrainCooling = 0;
	// For gas-grain collisional cooling, calculate the average dust temperature, weighted
	// by cross section. This will be the expected value of the temperatures a particle
	// 'experiences' when colliding with a grain.
	if (GASGRAINCOOL && numPop > 0)
	{
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
		// Cooling of the gas by collisions with the dust particles. This is a recipe
		// that does not depend on the grain present, except for the average temperature
		// (it probably assumes some average population, see coefficient in Krumholz et
		// al. (2011)). Only do this if there is dust, cause otherwise we dont have a
		// dust temperature of course.
		double Tdust = Tsum / weight;
		double Tgas = s.T;
		double nH = s.speciesNv[_inH];
		double nH2 = s.speciesNv[_inH2];
		gasGrainCooling = GasGrain::simpleGasGrainCool(Tdust, Tgas, nH, nH2);
	}
	return grainPhotoelectricHeating - gasGrainCooling;
}

double GasInterfaceImpl::lineCooling(const Solution& s) const
{
	// double result = _atomicLevels->cooling(s.HSolution);
	// if (_molecular)
	// 	result += _molecular->cooling(s.H2Solution);
	// return result;
	return 0.;
}

double GasInterfaceImpl::lineHeating(const Solution& s) const
{
	double result = _atomicLevels->netheating(s.HSolution);
	if (_molecular)
		result += _molecular->netheating(s.H2Solution);
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
		double dissheat = _molecular->dissociationHeating(s.H2Solution,
		                                                  s.specificIntensity);
		result += dissheat;
	}
	return result;
}
