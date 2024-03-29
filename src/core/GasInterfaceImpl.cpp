#include "GasInterfaceImpl.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "EigenAliases.hpp"
#include "Error.hpp"
#include "Functions.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "Ionization.hpp"
#include "Options.hpp"
#include "RadiationFieldTools.hpp"
#include "SimpleColumnFile.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "TwoPhoton.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <cmath>
#include <iostream>
#include <string>
#include <valarray>

namespace RADAGAST
{
    GasInterfaceImpl::GasInterfaceImpl(const Array& iFrequencyv, const Array& oFrequencyv, const Array& eFrequencyv)
        : _iFrequencyv{iFrequencyv}, _oFrequencyv{oFrequencyv}, _eFrequencyv{eFrequencyv}
    {
        // this file is simple enough to load and process it ad-hoc here
        InColumnFile h2cross_leiden("dat/h2/crossections_leiden.txt");
        h2cross_leiden.read(4, 1030);
        auto wav_nm = h2cross_leiden.column(0);
        // 1 absorption (ionization + LW); 2 LW (only dissociation, so ~10 times smaller than
        // absorption); 3 ionization only
        auto absorption_cs_cm2 = h2cross_leiden.column(3);

        int dataSize = wav_nm.size();
        Array frequencyv(dataSize);
        Array crossSectionv(dataSize);
        for (int i = 0; i < dataSize; i++)
        {
            int ri = dataSize - 1 - i;
            frequencyv[i] = Constant::LIGHT / (wav_nm[ri] * Constant::NM);
            crossSectionv[i] = absorption_cs_cm2[ri];
        }
        Spectrum data(frequencyv, crossSectionv);
        _h2crossv = data.binned(oFrequencyv);

        // Use approximation (constant * H cross section) for frequencies above the data range.
        // Otherwise, H2 will be transparent < 17 nm). Constant was determined by taking the
        // averate ratio with respect to the H cross section for all data < 28 nm.
        for (int i = oFrequencyv.size() - 1; i >= 0 && oFrequencyv[i] > data.freqMax(); i--)
            _h2crossv[i] = 3.33 * Ionization::crossSection(oFrequencyv[i]);
    }

    namespace
    {
        void sanitizeInput(double n, const std::valarray<double>& meanIntensityv, double fshield,
                           RADAGAST::GrainInterface& grainInfo)
        {
            Error::sanitize("n", n);
            Error::sanitize("fshield", fshield);
        }
    }

    void GasInterfaceImpl::updateGasState(RADAGAST::GasState& gs, double n, const std::valarray<double>& meanIntensityv,
                                          double fshield, RADAGAST::GrainInterface& grainInfo, GasDiagnostics* gd) const
    {
        sanitizeInput(n, meanIntensityv, fshield, grainInfo);
        Spectrum meanIntensity(_iFrequencyv, meanIntensityv);
        GasSolution s = solveTemperature(n, meanIntensity, fshield, grainInfo);
        s.setGasState(gs);
        if (gd) s.fillDiagnostics(gd);
    }

    Array GasInterfaceImpl::emissivityBasic(const RADAGAST::GasState& gs, bool SI) const
    {
        // There is some duplication from GasSolution::emissivityv, but let's keep both for now.
        Array emissivityv(_eFrequencyv.size());
        _freeBound.addEmissionCoefficientv(gs.temperature(), _eFrequencyv, emissivityv);
        _freeFree.addEmissionCoefficientv(gs.temperature(), _eFrequencyv, emissivityv);
        emissivityv *= npne(gs) / Constant::FPI;

        emissivityv += TwoPhoton::emissivityv(_eFrequencyv, gs._n2s);

        if (SI)
            return RadiationFieldTools::emissivity_to_SI(emissivityv);
        else
            return emissivityv;
    }

    Array GasInterfaceImpl::opacityBasic(const RADAGAST::GasState& gs, bool SI) const
    {
        auto sv = speciesVector(gs);

        Array opacityv(_oFrequencyv.size());
        _freeFree.addOpacityCoefficientv(gs.temperature(), _oFrequencyv, opacityv);
        opacityv *= sv.np() * sv.ne();

        for (size_t i = 0; i < _oFrequencyv.size(); i++)
            opacityv[i] += sv.nH() * Ionization::crossSection(_oFrequencyv[i]);

        // H2 ionization cross section is added separately
        opacityv += sv.nH2() * _h2crossv;

        // convert from cm-1 to m-1
        if (SI)
            return 100. * opacityv;
        else
            return opacityv;
    }

    Array GasInterfaceImpl::emissivityWithLines(const RADAGAST::GasState& gs, const Array& meanIntensityv,
                                                const GrainInterface& gri, bool SI, bool addHLines,
                                                bool addH2Lines) const
    {
        // get the continuum
        Array emissivityv = emissivityBasic(gs, false);

        // copy some things
        Spectrum meanIntensity(_iFrequencyv, meanIntensityv);
        auto sv = speciesVector(gs);

        // calculate the gas solution, without iteration since we already provide the correct t
        // and nv
        auto s = makeGasSolution(meanIntensity, 1., &gri);
        s.setT(gs._t);
        s.setSpeciesNv(sv.speciesNv());
        s.solveGrains();
        s.solveLevels();
        // TODO: technically we have to iterate here, since some H2 collision rates depend on
        // the ortho-para ratio, which is in turn derived from the levels populations. Ignoring
        // for now.

        if (addHLines) emissivityv += s.hModel()->emissivityv(_eFrequencyv);
        if (addH2Lines) emissivityv += s.h2Model()->emissivityv(_eFrequencyv);
        if (SI)
            return RadiationFieldTools::emissivity_to_SI(emissivityv);
        else
            return emissivityv;
    }

    Array GasInterfaceImpl::opacityWithLines(const RADAGAST::GasState& gs, const Array& meanIntensityv,
                                             const GrainInterface& gri, bool SI, bool addHLines, bool addH2Lines) const
    {
        // get the continuum
        Array opacityv = opacityBasic(gs, false);

        // copy some things
        Spectrum meanIntensity(_iFrequencyv, meanIntensityv);
        auto sv = speciesVector(gs);

        // calculate the gas solution, without iteration since we already provide the correct t
        // and nv
        auto s = makeGasSolution(meanIntensity, 1., &gri);
        s.setT(gs._t);
        s.setSpeciesNv(sv.speciesNv());

        // H2 formation on grains can influence H2 level populations
        s.solveGrains();
        s.solveLevels();

        // opacity of H lines
        if (addHLines) opacityv += s.hModel()->opacityv(_oFrequencyv);

        // opacity of H2 lines and continuum dissociation cross section (both zero when small H2
        // model is used)
        if (addH2Lines) opacityv += s.h2Model()->opacityv(_oFrequencyv);

        // convert from cm-1 to m-1
        if (SI)
            return 100. * opacityv;
        else
            return opacityv;
    }

    std::string GasInterfaceImpl::quickInfo(const RADAGAST::GasState& gs,
                                            const std::valarray<double>& meanIntensity) const
    {
        Spectrum si(_iFrequencyv, meanIntensity);
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
            const Spectrum* meanIntensity;
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
            double netHeat = p->gasInterfacePimpl->solveDensities(s, p->n, pow(10., logT), *p->meanIntensity,
                                                                  p->use_previous_solution);
            return netHeat;
        }
    }  // namespace

    GasSolution GasInterfaceImpl::solveTemperature(double n, const Spectrum& meanIntensity, double fshield,
                                                   RADAGAST::GrainInterface& gri) const
    {
        GasSolution s = makeGasSolution(meanIntensity, fshield, &gri);
        if (n <= 0)
        {
            solveDensities(s, 0, 0, meanIntensity);
            return s;
        }

        const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
        gsl_root_fsolver* solver = gsl_root_fsolver_alloc(T);

        gsl_function F;
        struct heating_f_params p = {this, n, &meanIntensity, &s, false};
        F.function = &heating_f;
        F.params = &p;

        const double Tmin = 1.;
        const double logTmin = log10(Tmin);
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
            std::cout << "Could not find equilibrium temperature. G0 is " << RadiationFieldTools::gHabing(meanIntensity)
                      << '\n';
            solveDensities(s, n, Tmaxmax, meanIntensity);
            return s;
        }
        if (heating_f_Tmin <= 0 && heating_f_Tmax <= 0)
        {
            // If even after all this, heating is negative on both sides, solve for the lowest
            // temperature as the 'safe' option.
            solveDensities(s, n, Tmin, meanIntensity);
            return s;
        }

        // Initialize the solver. Will fail if heating at Tmin and Tmax do not have different sign.
        int status = gsl_root_fsolver_set(solver, &F, logTmin, logTmax);
        if (status == GSL_EINVAL)
        {
            // If something else is stil wrong, be verbose and solve for lowest temperature
            std::cout << "heating at " << pow(10., logTmin) << " K -->" << heating_f_Tmin << '\n';
            std::cout << "heating at " << pow(10., logTmax) << " K -->" << heating_f_Tmax << '\n';
            solveDensities(s, n, Tmin, meanIntensity);
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

        double logTtolerance = std::log10(1. + Options::solvetemperature_Ttolerance);
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
            test_interval = gsl_root_test_interval(lower, upper, logTtolerance, 0);
            counter++;
        }
        gsl_root_fsolver_free(solver);
        // s was already updated during the last time GSL called the heating function
        return s;
    }

    GasSolution GasInterfaceImpl::solveDensities(double n, double T, const Spectrum& meanIntensity, double fshield,
                                                 RADAGAST::GrainInterface& gri, double h2FormationOverride) const
    {
        GasSolution s = makeGasSolution(meanIntensity, fshield, &gri);
        solveDensities(s, n, T, meanIntensity, false, h2FormationOverride);
        return s;
    }

    double GasInterfaceImpl::solveDensities(GasSolution& s, double n, double T, const Spectrum& meanIntensity,
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
            double ionFrac = Ionization::solveBalance(n, T, meanIntensity);
            double molFrac = 0.1;
            s.setSpeciesNv(guessSpeciesNv(n, ionFrac, molFrac));
        }
        // if nHtotal starts drifting, I should do someting extra here
        DEBUG("Set initial guess for speciesNv to\n" << s.speciesVector() << '\n');

        s.setT(T);

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
            // Update grain charge distributions and temperatures. Will affect the grain
            // photoelectric heating, the h2 formation rate, the h2 levels, and the gas-dust
            // energy exchange (gas cooling)
            s.solveGrains();

            // GRAIN AND CHEMISTRY SOLUTION -> SOURCE AND SINK RATES -> LEVEL SOLUTION
            s.solveLevels();

            if (previousAbundancev.array().isNaN().any())
            {
                std::cerr << previousAbundancev;
                Error::warn("Nan in chemistry solution! Setting to zero.");
                previousAbundancev.setZero();
                s.setSpeciesNv(EVector::Zero(previousAbundancev.size()));
            }

            // GRAIN, LEVELS, AND CHEMISTRY SOLUTIONS -> CHEM RATES
            double kFormH2 = 0;
            if (h2FormationOverride < 0)
                kFormH2 = s.kGrainH2FormationRateCoeff();
            else
                kFormH2 = h2FormationOverride;

            double kDissH2Levels = s.kDissH2Levels();
            if (kDissH2Levels < 0) Error::runtime("negative dissociation rate!");

            // Calculate all reaction rates and put them (+ the H2 rates provided as arguments)
            // into a vector of the right format
            EVector reactionRates = _chemistry.rateCoeffv(T, meanIntensity, kDissH2Levels, kFormH2);

            // CHEM RATES -> CHEMISTRY SOLUTION
            EVector newSpeciesNv = _chemistry.solveBalance(reactionRates, s.speciesVector().speciesNv());
            s.setSpeciesNv(newSpeciesNv);

            // TODO: add grain charging / photoelectric ejection to chemical network. Can be
            // important source of electrons sometimes.

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

    GasSolution GasInterfaceImpl::makeGasSolution(const Spectrum& meanIntensity, double fshield,
                                                  const RADAGAST::GrainInterface* gri) const
    {
        // Since I keep forgetting why it is good to make a factory function return unique
        // pointers, here is a reminder:

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
        std::unique_ptr<HModel> hm = _manager.makeHModel(&meanIntensity);
        std::unique_ptr<H2Model> h2m = _manager.makeH2Model(&meanIntensity, fshield);
        GasSolution s(gri, &meanIntensity, fshield, &_chemistry.speciesIndex(), move(hm), move(h2m), &_freeBound,
                      &_freeFree);
        return s;
    }

    EVector GasInterfaceImpl::guessSpeciesNv(double n, double ionToTotalFrac, double moleculeToNeutralFrac) const
    {
        return _chemistry.speciesIndex().linearCombination(SpeciesIndex::e_p_H_H2,
                                                           {ionToTotalFrac * n, ionToTotalFrac * n,
                                                            (1 - moleculeToNeutralFrac) * (1 - ionToTotalFrac) * n,
                                                            moleculeToNeutralFrac * (1 - ionToTotalFrac) * n / 2});
    }

    double GasInterfaceImpl::npne(const RADAGAST::GasState& gs) const
    {
        return gs.density(_chemistry.speciesIndex().index("H+")) * gs.density(_chemistry.speciesIndex().index("e-"));
    }

    double GasInterfaceImpl::nH(const RADAGAST::GasState& gs) const
    {
        return gs.density(_chemistry.speciesIndex().index("H"));
    }

    SpeciesVector GasInterfaceImpl::speciesVector(const RADAGAST::GasState& gs) const
    {
        SpeciesVector sv(&_chemistry.speciesIndex());
        sv.setDensities(Eigen::Map<const EVector>(&gs._nv[0], gs._nv.size()));
        return sv;
    }
}
