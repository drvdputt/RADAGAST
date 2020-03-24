#include "GasInterfaceImpl.hpp"
// #include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "RadiationFieldTools.hpp"
#include "SimpleColumnFile.hpp"
#include "SpecialFunctions.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "TwoPhoton.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <iostream>

namespace GasModule
{
    GasInterfaceImpl::GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv, const Array& eFrequencyv,
                                       const std::string& atomChoice, const std::string& moleculeChoice)
        : _iFrequencyv{iFrequencyv}, _oFrequencyv{oFrequencyv}, _eFrequencyv{eFrequencyv},
          _manager(atomChoice, moleculeChoice)
    {
        InColumnFile h2cross_leiden("dat/h2/crossections_leiden.txt");
        h2cross_leiden.read(4, 1030);
        auto wav_nm = h2cross_leiden.column(0);
        auto absorption_cs_cm2 = h2cross_leiden.column(1);

        int dataSize = wav_nm.size();
        Array frequencyv(dataSize);
        Array crossSectionv(dataSize);
        for (int i = 0; i < dataSize; i++)
        {
            int ri = dataSize - 1 - i;
            frequencyv[i] = Constant::LIGHT / (wav_nm[ri] * Constant::NM);
            crossSectionv[i] = absorption_cs_cm2[ri];
        }
        _h2crossv = Spectrum(frequencyv, crossSectionv).binned(oFrequencyv);
    }

    void GasInterfaceImpl::updateGasState(GasModule::GasState& gs, double n,
                                          const std::valarray<double>& specificIntensityv,
                                          GasModule::GrainInterface& grainInfo, GasDiagnostics* gd) const
    {
        Spectrum specificIntensity(_iFrequencyv, specificIntensityv);
        GasSolution s = solveTemperature(n, specificIntensity, grainInfo);
        s.setGasState(gs);
        if (gd) s.fillDiagnostics(gd);
    }

    Array GasInterfaceImpl::emissivity(const GasModule::GasState& gs, bool SI) const
    {
        // There is some duplication from GasSolution::emissivityv, but let's keep both for now.
        Array emissivityv(_eFrequencyv.size());
        _freeBound.addEmissionCoefficientv(gs.temperature(), _eFrequencyv, emissivityv);
        _freeFree.addEmissionCoefficientv(gs.temperature(), _eFrequencyv, emissivityv);
        emissivityv *= npne(gs) / Constant::FPI;

        // TODO: switch to add line emissivity (needs recalculation of levels)

        emissivityv += TwoPhoton::emissivityv(_eFrequencyv, gs._n2s);

        if (SI)
            return RadiationFieldTools::emissivity_to_SI(emissivityv);
        else
            return emissivityv;
    }

    Array GasInterfaceImpl::opacity(const GasModule::GasState& gs, bool SI) const
    {
        auto sv = speciesVector(gs);

        Array opacityv(_oFrequencyv.size());
        _freeFree.addOpacityCoefficientv(gs.temperature(), _oFrequencyv, opacityv);
        opacityv *= sv.np() * sv.ne();

        // TODO: this should actually be the average over the cross section for this frequency bin
        for (size_t i = 0; i < _oFrequencyv.size(); i++)
            opacityv[i] += sv.nH() * Ionization::crossSection(_oFrequencyv[i]);

        bool withLines = false;
        if (withLines)
        {
            // TODO: make H model, calculate balance and get lines (no iteration needed if the
            // given gas state is already self-consistent).

            // TODO: H2 opacity. Maybe use smoothed data file from Heays et al. 2017 (data on E. van
            // Dishoeck's home page)

            // need grains here for h2 formation, which in turn I need for h2model.solve()

            // auto h2model = _manager.makeH2Model();
            // auto sv = speciesVector(gs);
            // CollisionParameters cp(gs._t, sv);
            // h2model.solve(sv.nH2(), cp, specificIntensity, h2form);
        }

        if (!withLines) opacityv += sv.nH2() * _h2crossv;

        // convert from cm-1 to m-1
        if (SI)
            return 100. * opacityv;
        else
            return opacityv;
    }

    std::string GasInterfaceImpl::quickInfo(const GasModule::GasState& gs,
                                            const std::valarray<double>& specificIntensity) const
    {
        Spectrum si(_iFrequencyv, specificIntensity);
        auto sv = speciesVector(gs);
        std::stringstream ss;
        ss << "G0 " << RadiationFieldTools::gHabing(si) << " T " << gs._t << " ne " << sv.ne() << " np " << sv.np()
           << " nH " << sv.nH() << " nH2 " << sv.nH2();
        return ss.str();
    }

    int GasInterfaceImpl::index(const std::string& name) const { return _chemistry.speciesIndex().index(name); }

    namespace
    {
        struct heating_f_params
        {
            const GasInterfaceImpl* gasInterfacePimpl;
            double n;
            const Spectrum* specificIntensity;
            GasSolution* solution_storage;
            bool use_previous_solution;
        };

        // Function that will be used by the GSL search algorithm to find the equilibrium temperature.
        // The state of the system will be updated every time the algorithm calls this function (and
        // stored via the pointer suppied in the heating_f_params struct).
        double heating_f(double logT, void* params)
        {
            auto* p = static_cast<struct heating_f_params*>(params);
            // Refresh the solution stored somewhere, with optimization based on the current solution
            GasSolution& s = *p->solution_storage;
            double netHeat = p->gasInterfacePimpl->solveDensities(s, p->n, pow(10., logT), *p->specificIntensity,
                                                                  p->use_previous_solution);
            return netHeat;
        }
    }  // namespace

    GasSolution GasInterfaceImpl::solveTemperature(double n, const Spectrum& specificIntensity,
                                                   GasModule::GrainInterface& gri) const
    {
        GasSolution s = makeGasSolution(specificIntensity, &gri);
        if (n <= 0)
        {
            solveDensities(s, 0, 0, specificIntensity);
            return s;
        }

        const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
        gsl_root_fsolver* solver = gsl_root_fsolver_alloc(T);

        gsl_function F;
        struct heating_f_params p = {this, n, &specificIntensity, &s, false};
        F.function = &heating_f;
        F.params = &p;

        const double Tmin = 1.;
        const double logTmin = log10(Tmin);
        const double logTtolerance = 1.e-3;
        // Assume that the heating is positive at Tmin (this has yet to fail). Then, find a suitable
        // upper limit for the bracket (one where the heating is negative. Try 10000 first, which
        // should be suitable for most applications. Then, keep multiplying the temperature with a
        // constant factor to see if it helps.
        const double Tmax = 10000.;
        const double Tmaxmax = 1.e7;
        double logTmax = log10(Tmax);
        const double logTmaxmax = log10(Tmaxmax);
        double heating_f_Tmax = heating_f(logTmax, &p);
        double heating_f_Tmin = heating_f(logTmin, &p);
        double logFactor = log10(3.);
        while (heating_f_Tmax >= 0. && logTmax < logTmaxmax)
        {
            logTmax += logFactor;
            heating_f_Tmax = heating_f(logTmax, &p);
        }
        if (logTmax > logTmaxmax && heating_f_Tmax >= 0.)
        {
            // If net heating is still positive, just solve for the maximum temperature.
            std::cout << "Could not find equilibrium temperature. G0 is "
                      << RadiationFieldTools::gHabing(specificIntensity) << '\n';
            solveDensities(s, n, Tmaxmax, specificIntensity);
            return s;
        }

        // Initialize the solver. Will fail if heating at Tmin and Tmax do not have different sign.
        int status = gsl_root_fsolver_set(solver, &F, logTmin, logTmax);
        if (status == GSL_EINVAL)
        {
            // The heating is most likely negative on both sides in this case. Just solve for the
            // lowest temperature as the 'safe' option.
            std::cout << "heating at " << pow(10., logTmin) << " K -->" << heating_f_Tmin << '\n';
            std::cout << "heating at " << pow(10., logTmax) << " K -->" << heating_f_Tmax << '\n';
            solveDensities(s, n, Tmin, specificIntensity);
            return s;
        }

        // Iterate once to initialize the solution
        status = gsl_root_fsolver_iterate(solver);
        if (status)
        {
            std::cerr << gsl_strerror(status) << " in gsl root iteration\n";
            return s;
        }
        // Then, start using the current solution as an initial guess of the next one
        p.use_previous_solution = true;

        int test_interval = GSL_CONTINUE;
        int counter = 0;
        while (test_interval != GSL_SUCCESS)
        {
            status = gsl_root_fsolver_iterate(solver);
            if (status)
            {
                std::cerr << gsl_strerror(status) << " in gsl root iteration\n";
                return s;
            }
            double lower = gsl_root_fsolver_x_lower(solver);
            double upper = gsl_root_fsolver_x_upper(solver);
            test_interval = gsl_root_test_interval(lower, upper, logTtolerance, logTtolerance);
            counter++;
        }
        gsl_root_fsolver_free(solver);
        // s was already updated during the last time GSL called the heating function
        return s;
    }

    GasSolution GasInterfaceImpl::solveDensities(double n, double T, const Spectrum& specificIntensity,
                                                 GasModule::GrainInterface& gri, double h2FormationOverride) const
    {
        GasSolution s = makeGasSolution(specificIntensity, &gri);
        solveDensities(s, n, T, specificIntensity, false, h2FormationOverride);
        return s;
    }

    double GasInterfaceImpl::solveDensities(GasSolution& s, double n, double T, const Spectrum& specificIntensity,
                                            bool startFromCurrent, double h2FormationOverride) const
    {
        if (n <= 0)
        {
            s.setT(T);
            s.makeZero();
            return 0.;
        }

        DEBUG("Calculating densities for T = " << T << "K\n");

        // Decide how to do the initial guess; manual, or using the previous state.
        bool manualGuess = true;

        // If a previous solution is available (pointer is nonzero), then check if we want to use it
        if (startFromCurrent)
        {
            double fT = T / s.t();
            if (fT > 0.75 && fT < 1.5)
            {
                DEBUG("Using previous speciesNv as initial guess\n");
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
        // if nHtotal starts drifting, I should do someting extra here
        DEBUG("Set initial guess for speciesNv to\n" << s.speciesVector() << '\n');

        s.setT(T);

        // Formation rate on grain surfaces. We declare it here so that it can be passed to the level
        // population calculation too (needed for H2 formation pumping).
        double kFormH2;
        if (h2FormationOverride >= 0)
            kFormH2 = h2FormationOverride;
        else
            kFormH2 = s.kGrainH2FormationRateCoeff();

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
        EVector previousAbundancev = s.speciesVector().speciesNv();
        double previousHeating = 0;
        double previousCooling = 0;
        while (!stopCriterion)
        {
            // CHEMISTRY SOLUTION -> SOURCE AND SINK RATES -> LEVEL SOLUTION
            s.solveLevels(kFormH2 * s.speciesVector().nH());

            if (previousAbundancev.array().isNaN().any())
            {
                std::cerr << previousAbundancev;
                Error::runtime("Nan in chemistry solution!");
            }

            // LEVELS AND CHEMISTRY SOLUTIONS -> CHEM RATES
            if (h2FormationOverride < 0)
            {
                // update h2 formation rate if not overriden by a constant (need to update every
                // iteration because grain temperature can change)
                kFormH2 = s.kGrainH2FormationRateCoeff();
            }

            double kDissH2Levels = s.kDissH2Levels();
            if (kDissH2Levels < 0) Error::runtime("negative dissociation rate!");

            // CHEM RATES -> CHEMISTRY SOLUTION
            EVector reactionRates = _chemistry.rateCoeffv(T, specificIntensity, kDissH2Levels, kFormH2);
            EVector newSpeciesNv = _chemistry.solveBalance(reactionRates, s.speciesVector().speciesNv());
            s.setSpeciesNv(newSpeciesNv);

            // Recalculate the grain charge distributions and temperatures, which will affect the grain
            // photoelectric heating, the h2 formation rate, and the gas-dust energy exchange (gas
            // cooling)
            s.solveGrains();
            // TODO: Add effect of grain charging to chemical network. I think it might be possible to
            // do this by imposing a conservation equation for the number of electrons: ne + nH + nH2 =
            // (ne + nH + nH2)_0 + <Cg>*ng. Another option would be to include the grain charge rates
            // into the network as extra reactions. Grain recombination / charge exchange reactions
            // could also be added.

            // CONVERGENCE CHECK. An abundance has converged if it changes by less than 1%, or if the
            // total amount is negligible compared to the norm.
            EArray changev = newSpeciesNv - previousAbundancev;
            double norm = newSpeciesNv.norm();
            Eigen::Array<bool, Eigen::Dynamic, 1> convergedv =
                changev.abs() <= 1e-3 * previousAbundancev.array() || newSpeciesNv.array() < 1.e-30 * norm;

            double newHeating = s.heating();
            double newCooling = s.cooling();
            bool heatingConverged = TemplatedUtils::equalWithinTolerance(newHeating, previousHeating, 1e-2);
            bool coolingConverged = TemplatedUtils::equalWithinTolerance(newCooling, previousCooling, 1e-2);
            counter++;
            DEBUG("Chemistry: " << counter << '\n'
                                << s.speciesVector() << '\n'
                                << "convergence: \n"
                                << convergedv << '\n');
            DEBUG("New heat: " << newHeating << " previous: " << previousHeating << '\n');
            DEBUG("New cool: " << newCooling << " previous: " << previousCooling << '\n');

            // Decide if we want to stop
            bool allQuantitiesConverged = convergedv.all() && heatingConverged && coolingConverged;
            stopCriterion = allQuantitiesConverged || counter > Options::densities_maxiterations;

            // Save these results so they can be used for the next convergence check
            previousAbundancev = newSpeciesNv;
            previousHeating = newHeating;
            previousCooling = newCooling;
        }
        return previousHeating - previousCooling;
    }

    GasSolution GasInterfaceImpl::makeGasSolution(const Spectrum& specificIntensity,
                                                  const GasModule::GrainInterface* gri) const
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
        GasSolution s(gri, specificIntensity, &_chemistry.speciesIndex(), move(hm), move(h2m), _freeBound, _freeFree);
        return s;
    }

    EVector GasInterfaceImpl::guessSpeciesNv(double n, double ionToTotalFrac, double moleculeToNeutralFrac) const
    {
        return _chemistry.speciesIndex().linearCombination(SpeciesIndex::e_p_H_H2,
                                                           {ionToTotalFrac * n, ionToTotalFrac * n,
                                                            (1 - moleculeToNeutralFrac) * (1 - ionToTotalFrac) * n,
                                                            moleculeToNeutralFrac * (1 - ionToTotalFrac) * n / 2});
    }

    double GasInterfaceImpl::npne(const GasModule::GasState& gs) const
    {
        return gs.density(_chemistry.speciesIndex().index("H+")) * gs.density(_chemistry.speciesIndex().index("e-"));
    }

    double GasInterfaceImpl::nH(const GasModule::GasState& gs) const
    {
        return gs.density(_chemistry.speciesIndex().index("H"));
    }

    SpeciesVector GasInterfaceImpl::speciesVector(const GasModule::GasState& gs) const
    {
        SpeciesVector sv(&_chemistry.speciesIndex());
        sv.setDensities(Eigen::Map<const EVector>(&gs._nv[0], gs._nv.size()));
        return sv;
    }
}
