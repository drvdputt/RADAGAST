#include "GasInterfaceImpl.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GasGrainInteraction.hpp"
#include "GrainPhotoelectricEffect.hpp"
#include "GrainType.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "Timer.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

using namespace std;

GasInterfaceImpl::GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv,
                                   const Array& eFrequencyv, const string& atomChoice,
                                   const string& moleculeChoice)
                : _iFrequencyv{iFrequencyv}, _oFrequencyv{oFrequencyv},
                  _eFrequencyv{eFrequencyv}, _manager(atomChoice, moleculeChoice)
{
}

void GasInterfaceImpl::updateGasState(GasModule::GasState& gs, double n,
                                      const valarray<double>& specificIntensityv,
                                      GasModule::GrainInterface& grainInfo,
                                      GasDiagnostics* gd) const
{
	Timer t("Update gas state");
	Spectrum specificIntensity(_iFrequencyv, specificIntensityv);
	GasSolution s = solveTemperature(n, specificIntensity, grainInfo);
	gs = s.makeGasState(_oFrequencyv, _eFrequencyv);
	if (gd)
		s.fillDiagnostics(gd);
}

void GasInterfaceImpl::initializeGasState(GasModule::GasState& gs, double n, double T,
                                          GasModule::GrainInterface& gri,
                                          GasDiagnostics* gd) const
{
	GasSolution s = solveInitialGuess(n, T, gri);
	gs = s.makeGasState(_oFrequencyv, _eFrequencyv);
	if (gd)
		s.fillDiagnostics(gd);
}

double GasInterfaceImpl::emissivity_SI(const GasModule::GasState& gs, size_t iFreq) const
{
	return 0.1 * gs._emissivityv[iFreq];
}

double GasInterfaceImpl::opacity_SI(const GasModule::GasState& gs, size_t iFreq) const
{
	return 100 * gs._opacityv[iFreq];
}

GasSolution GasInterfaceImpl::solveInitialGuess(double n, double T,
                                                GasModule::GrainInterface& gri) const
{
	Array iFrequencyv = Testing::defaultFrequencyv(1000);
	Array blackbodyv(iFrequencyv.size());
	for (size_t i = 0; i < iFrequencyv.size(); i++)
		blackbodyv[i] = SpecialFunctions::planck(iFrequencyv[i], T);
	return solveTemperature(n, Spectrum(iFrequencyv, blackbodyv), gri);
}

namespace
{
struct heating_f_params
{
	const GasInterfaceImpl* gasInterfacePimpl;
	double n;
	const Spectrum* specificIntensity;
	GasModule::GrainInterface* grainInterface;
	GasSolution* solution_storage;
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
	GasSolution& s = *p->solution_storage;
	p->gasInterfacePimpl->solveDensities(s, p->n, pow(10., logT), *p->specificIntensity,
	                                     *p->grainInterface, p->use_previous_solution);
	double heat = s.heating();
	double cool = s.cooling();
	return heat - cool;
}
} // namespace

GasSolution GasInterfaceImpl::solveTemperature(double n, const Spectrum& specificIntensity,
                                               GasModule::GrainInterface& gri) const
{
	GasSolution s = makeGasSolution(specificIntensity, gri);
	if (n <= 0)
		solveDensities(s, 0, 0, specificIntensity, gri);
	else
	{
		const double Tmax = 100000.;
		const double Tmin = 1.;
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

GasSolution GasInterfaceImpl::solveDensities(double n, double T,
                                             const Spectrum& specificIntensity,
                                             GasModule::GrainInterface& gri,
                                             double h2FormationOverride) const
{
	GasSolution s = makeGasSolution(specificIntensity, gri);
	solveDensities(s, n, T, specificIntensity, gri, false, h2FormationOverride);
	return s;
}

void GasInterfaceImpl::solveDensities(GasSolution& s, double n, double T,
                                      const Spectrum& specificIntensity,
                                      GasModule::GrainInterface& gri, bool startFromCurrent,
                                      double h2FormationOverride) const
{
	if (n <= 0)
	{
		s.setT(T);
		s.makeZero();
		return;
	}

	DEBUG("Calculating densities for T = " << T << "K" << endl);

	// Decide how to do the initial guess; manual, or using the previous state.
	bool manualGuess = true;

	// If a previous solution is available (pointer is nonzero), then check if we want to
	// use it
	if (startFromCurrent)
	{
		double fT = T / s.t();
		if (fT > 0.75 && fT < 1.5)
		{
			DEBUG("Using previous speciesNv as initial guess" << std::endl);
			manualGuess = false;
		}
		// else manualGuess stays true
	}
	if (manualGuess)
	{
		double ionFrac = Ionization::solveBalance(n, T, specificIntensity);
		double molFrac = 0.1;
		s.setSpeciesNv(guessSpeciesNv(n, ionFrac, molFrac));
	}

	DEBUG("Set initial guess for speciesNv to\n" << s.speciesNv() << '\n');

	// else
	// {
	// 	// make sure that nHtotal and the initial guess are consistent
	// 	double nNeutralOfSpeciesNv = s.speciesNv()(_inH) + 2 * s.speciesNv()(_inH2);
	// 	double nHtotalOfSpeciesNv = s.speciesNv()(_inp) + nNeutralOfSpeciesNv;
	// 	double ionFrac = s.speciesNv()(_inp) / nHtotalOfSpeciesNv;
	// 	double molFrac = 2 * s.speciesNv()(_inH2) / nNeutralOfSpeciesNv;
	// 	s.setSpeciesNv(guessSpeciesNv(n, ionFrac, molFrac));
	// }
	s.setT(T);

	// Formation rate on grain surfaces. We declare it here so that it can be passed to the
	// level population calculation too (needed for H2 formation pumping).
	double kFormH2;
	if (h2FormationOverride >= 0)
		kFormH2 = h2FormationOverride;
	else
		kFormH2 = GasGrain::surfaceH2FormationRateCoeff(gri, T);

	/* The main loop:
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
	EVector previousAbundancev = s.speciesNv();
	double previousHeating = 0;
	double previousCooling = 0;
	while (!stopCriterion)
	{
		// CHEMISTRY SOLUTION -> SOURCE AND SINK RATES -> LEVEL SOLUTION
		s.solveLevels(kFormH2 * s.speciesNv()(_inH));

		if (previousAbundancev.array().isNaN().any())
			Error::runtime("Nan in chemistry solution!");

		// LEVELS AND CHEMISTRY SOLUTIONS -> CHEM RATES
		if (h2FormationOverride < 0)
		{
			// update h2 formation rate if not overriden by a constant (need to
			// update every iteration because grain temperature can change)
			kFormH2 = GasGrain::surfaceH2FormationRateCoeff(gri, T);
		}

		double kDissH2Levels = s.kDissH2Levels();
		if (kDissH2Levels < 0)
			Error::runtime("negative dissociation rate!");

		// CHEM RATES -> CHEMISTRY SOLUTION
		EVector reactionRates = _chemistry.rateCoeffv(T, specificIntensity, kDissH2Levels, kFormH2);
		EVector newSpeciesNv = _chemistry.solveBalance(reactionRates, s.speciesNv());
		s.setSpeciesNv(newSpeciesNv);

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

		// CONVERGENCE CHECK. An abundance has converged if it changes by less than 1%,
		// or if the total amount is negligible compared to the norm.
		EArray changev = s.speciesNv() - previousAbundancev;
		double norm = s.speciesNv().norm();
		Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
		                changev.abs() <= 1e-3 * previousAbundancev.array() ||
		                s.speciesNv().array() < 1.e-99 * norm;

		// Changing the grain temperature influences the following:
		// h2 formation rate
		// gas-dust energy exchange (gas cooling)
		updateGrainTemps(s, gri);

		double newHeating = s.heating();
		double newCooling = s.cooling();
		bool heatingConverged = TemplatedUtils::equalWithinTolerance(
		                newHeating, previousHeating, 1e-2);
		bool coolingConverged = TemplatedUtils::equalWithinTolerance(
		                newCooling, previousCooling, 1e-2);
		counter++;
		DEBUG("Chemistry: " << counter << '\n'
		                    << s.speciesNv() << '\n'
		                    << "convergence: \n"
		                    << convergedv << '\n');
		DEBUG("New heat: " << newHeating << " previous: " << previousHeating << '\n');
		DEBUG("New cool: " << newCooling << " previous: " << previousCooling << '\n');

		previousAbundancev = s.speciesNv();
		previousHeating = newHeating;
		previousCooling = newCooling;

		bool allQuantitiesConverged =
		                convergedv.all() && heatingConverged && coolingConverged;

		// Currently, the implementation without molecules does not need iteration.
		stopCriterion = allQuantitiesConverged ||
		                counter > Options::densities_maxiterations;
	}
}

// GasSolution GasInterfaceImpl::solveDensitiesNoH2(double n, double T,
//                                                  const Spectrum& specificIntensity,
//                                                  const GasModule::GrainInterface&) const;
// {
// // When ignoring H2
// // DIRECT SOLUTION (IONIZATION BALANCE ONLY, NO MOLECULES)

// // Just solve the ionization balance in the nebular approximation.
// double f = Ionization::solveBalance(nHtotal, T, specificIntensity);
// DEBUG("Ionized fraction = " << f << endl);

// // Neutral fraction
// s.speciesNv[_inH] = nHtotal * (1 - f);
// // Ionized fraction
// s.speciesNv[_inp] = nHtotal * f;
// // Electron density is simply equal to proton density
// s.speciesNv[_ine] = s.speciesNv[_inp];
// s.speciesNv[_inH2] = 0;
// }

void GasInterfaceImpl::updateGrainTemps(const GasSolution& s,
                                        GasModule::GrainInterface& g) const
{
	GrainPhotoelectricEffect::Environment env(
	                s.specificIntensity(), s.t(), s.ne(), s.np(), {-1, 1, 0, 0},
	                {s.ne(), s.np(), s.nH(), s.nH2()},
	                {Constant::ELECTRONMASS, Constant::PROTONMASS, Constant::HMASS_CGS,
	                 2 * Constant::HMASS_CGS});

	for (auto& pop : *g.populationv())
	{
		const GrainType* type = pop.type();

		size_t numSizes = pop.numSizes();
		Array grainHeatPerSizev(numSizes);
		// Array grainPhotoPerSizev(numSizes);

		if (type->heatingAvailable() && Options::cooling_gasGrainCollisions)
		{
			GrainPhotoelectricEffect gpe(*type);

			for (int m = 0; m < numSizes; m++)
			{
				auto cd = gpe.calculateChargeDistribution(pop.size(m), env,
				                                          pop.qAbsv(m));
				grainHeatPerSizev[m] = gpe.gasGrainCollisionCooling(
				                pop.size(m), env, cd, pop.temperature(m), true);
				// cout << "extra grain heat " << m << " " << grainHeatPerSizev[m] << '\n';

				// grainPhotoPerSizev[m] = gpe.heatingRateA(pop.size(m), env, pop.qAbsv(m), cd);
				// cout << "- photo heat " << m << " " << grainPhotoPerSizev[m] << '\n';
			}
		}

		const Array& h2FormationHeatv =
		                GasGrain::surfaceH2FormationHeatPerSize(pop, s.t(), s.nH());

		pop.recalculateTemperature(s.specificIntensity().frequencyv(),
		                           s.specificIntensity().valuev(),
		                           grainHeatPerSizev + h2FormationHeatv);
	}
}

GasSolution GasInterfaceImpl::makeGasSolution(const Spectrum& specificIntensity,
                                              const GasModule::GrainInterface& gri) const
{
	// TODO: rethink how the GasSolution is passed around, taking the following into account:

	// std::unique_ptr is the C++11 way to express exclusive ownership, but one of its most
	// attractive features is that it easily and efficiently converts to a std::shared_ptr.

	// This is a key part of why std::unique_ptr is so well suited as a factory function
	// return type. Factory functions can’t know whether callers will want to use exclusive
	// ownership semantics for the object they return or whether shared ownership (i.e.,
	// std::shared_ptr) would be more appropriate. By returning a std::unique_ptr, factories
	// provide callers with the most efficient smart pointer, but they don’t hinder callers
	// from replacing it with its more flexible sibling.

	// std::shared_ptr to std::unique_ptr is not allowed. Once you’ve turned lifetime
	// management of a resource over to a std::shared_ptr, there’s no changing your mind.
	// Even if the reference count is one, you can’t reclaim ownership of the resource in
	// order to, say, have a std::unique_ptr manage it.

	// Reference: Effective Modern C++. 42 SPECIFIC WAYS TO IMPROVE YOUR USE OF C++11 AND
	// C++14. Scott Meyers.
	std::unique_ptr<HModel> hm = _manager.makeHModel();
	std::unique_ptr<H2Model> h2m = _manager.makeH2Model();
	GasSolution s(gri, specificIntensity, move(hm), move(h2m), _freeBound, _freeFree);
	return s;
}

EVector GasInterfaceImpl::guessSpeciesNv(double n, double ionToTotalFrac,
                                         double moleculeToNeutralFrac) const
{
	EVector speciesNv = EVector::Zero((SpeciesIndex::size()));
	speciesNv(_ine) = ionToTotalFrac * n;
	speciesNv(_inp) = speciesNv(_ine);
	speciesNv(_inH) = (1 - moleculeToNeutralFrac) * (1 - ionToTotalFrac) * n;
	speciesNv(_inH2) = moleculeToNeutralFrac * (1 - ionToTotalFrac) * n / 2;
	return speciesNv;
}
