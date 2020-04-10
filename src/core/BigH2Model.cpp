#include "BigH2Model.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "H2Data.hpp"
#include "LevelSolver.hpp"
#include <sstream>

namespace GasModule
{
    void BigH2Model::solve(double n, const CollisionParameters& cp, const Spectrum& specificIntensity, double h2form)
    {
        _levelSolution.setT(cp._t);
        _n = n;

        if (n <= 0)
        {
            int numLv = _h2Data->numLv();
            _levelSolution.setCvv(EMatrix::Zero(numLv, numLv));
            _levelSolution.setNv(EVector::Zero(numLv));
            return;
        }

        _tvv = _h2Data->totalTransitionRatesvv(specificIntensity, cp, &_cvv, &_bpvv);
        _levelSolution.setCvv(_cvv);

        // TODO: Use better formation pumping recipe
        EVector sourcev = EVector::Zero(_h2Data->numLv());
        sourcev.head(_h2Data->startOfExcitedIndices()) = _h2Data->formationDistribution();
        sourcev *= h2form / sourcev.sum();

        EVector sinkv = dissociationSinkv(specificIntensity);

        if (!_levelSolution.isNvSet())
        {
            EVector initialGuessv = EVector::Zero(_h2Data->numLv());

            // Use LTE for the X levels, and 0 for the rest
            int endX = _h2Data->startOfExcitedIndices();
            initialGuessv.head(endX) = LevelSolver::statisticalEquilibrium_boltzman(n, cp._t, _h2Data->ev().head(endX),
                                                                                    _h2Data->gv().head(endX));

            DEBUG("Using LTE as initial guess for H2" << std::endl);
            _levelSolution.setNv(initialGuessv);
        }

        int fullyConnectedCutoff = _h2Data->startOfExcitedIndices();
        EVector newNv = LevelSolver::statisticalEquilibrium_iterative(n, _tvv, sourcev, sinkv, _levelSolution.nv(),
                                                                      fullyConnectedCutoff);
        _levelSolution.setNv(newNv);
    }

    double BigH2Model::dissociationRate(const Spectrum& specificIntensity) const
    {
#ifdef STERNBERG2014
        // TODO: move this to simple h2 model
        // See 2014-Sternberg eq 3
        auto iv = specificIntensity.valuev();
        auto nuv = specificIntensity.frequencyv();

        // F0 = integral 912 to 1108 Angstrom of Fnu(= 4pi Inu) with Inu in cm-2 s-1 Hz sr-1
        Array photonFluxv = Constant::FPI * iv / nuv / Constant::PLANCK;
        constexpr double freqLWmin{Constant::LIGHT / 1108 / Constant::ANGSTROM};
        constexpr double freqLWmax{Constant::LIGHT / 912 / Constant::ANGSTROM};
        size_t iLWmin{TemplatedUtils::index(freqLWmin, nuv)};
        size_t iLWmax{TemplatedUtils::index(freqLWmax, nuv)};
        double F0 = TemplatedUtils::integrate<double>(nuv, photonFluxv, iLWmin, iLWmax);

        // eq 4 and 5
        double Iuv{F0 / 2.07e7};
        double result{5.8e-11 * Iuv};
        return result;
#else
        EVector directv = directDissociationIntegralv(specificIntensity);
        EVector solomonv = spontaneousDissociationSinkv();
        EVector fv = _levelSolution.fv();

        // Dot product = total rate [cm-3 s-1]. Divide by total to get [s-1] rate, which
        // can be used in chemical network (it will multiply by the density again).
        double directFractional = directv.dot(fv);
        double solomonFractional = solomonv.dot(fv);
        DEBUG("Dissociation: direct rate:" << directFractional << " solomon rate: " << solomonFractional << '\n');
        return directFractional + solomonFractional;
#endif
    }

    double BigH2Model::dissociationHeating(const Spectrum& specificIntensity) const
    {
        // Fraction that dissociates per second * kinetic energy per dissociation * density of
        // level population = heating power density
        auto p = _h2Data->dissociationProbabilityv().array();
        auto k = _h2Data->dissociationKineticEnergyv().array();
        double solomonHeat = _levelSolution.nv().dot((p * k).matrix());

        // TODO: precalculate this?
        EVector directHeatv = directDissociationIntegralv(specificIntensity, true);
        double directHeat = _levelSolution.nv().dot(directHeatv);

        DEBUG("Dissociation heat: direct heat: " << directHeat << " solomon heat:" << solomonHeat << '\n');
        return solomonHeat + directHeat;
    }

    double BigH2Model::netHeating() const
    {
        // TODO cooling effect of collisional dissociation?
        return _levelSolution.netHeating();
    }

    double BigH2Model::orthoPara() const
    {
        if (!_levelSolution.isNvSet()) return 0.75;

        double orthoSum = 0;
        double paraSum = 0;
        for (size_t i = 0; i < _h2Data->numLv(); i++)
        {
            double n = _levelSolution.nv()(i);
            if (_h2Data->level(i).ortho())
                orthoSum += n;
            else
                paraSum += n;
        }
        return orthoSum / (orthoSum + paraSum);
    }

    Array BigH2Model::emissivityv(const Array& eFrequencyv) const { return _levelSolution.emissivityv(eFrequencyv); }

    Array BigH2Model::opacityv(const Array& oFrequencyv) const
    {
        // Start with the line opacity
        Array totalOpv = _levelSolution.opacityv(oFrequencyv);

        // Then add the dissociation cross sections of each level
        for (size_t iLv : _h2Data->levelsWithCrossSectionv())
        {
            const std::vector<Spectrum>& csv = _h2Data->directDissociationCrossSections(iLv);
            for (const Spectrum& cs : csv) totalOpv += cs.binned(oFrequencyv);
            // TODO: Discovered bug. Need to multiply with level density here!
        }
        return totalOpv;
    }

    void BigH2Model::extraDiagnostics(GasDiagnostics& gd, const Spectrum& specificIntensity) const
    {
        auto levelLabel = [](const H2Data::H2Level& level) {
            std::stringstream label;
            label << "E" << static_cast<int>(level.eState()) << " J" << level.j() << " v" << level.v();
            return label.str();
        };

        EVector fv = _levelSolution.fv();

        EVector contDissv = directDissociationIntegralv(specificIntensity);
        gd.setUserValue("H2 contdiss", contDissv.dot(fv));
        for (int i : _h2Data->levelsWithCrossSectionv())
            gd.setUserValue("H2 contdiss " + levelLabel(_h2Data->level(i)), contDissv[i] * fv[i]);

        EVector solomonDissv = spontaneousDissociationSinkv();
        gd.setUserValue("H2 solomon", solomonDissv.dot(fv));
        for (int i = 0; i < _h2Data->numLv(); i++)
        {
            double dissContribution = solomonDissv[i] * fv[i];
            if (dissContribution > 0) gd.setUserValue("H2 solomon " + levelLabel(_h2Data->level(i)), dissContribution);
        }
    }

    EVector BigH2Model::dissociationSinkv(const Spectrum& specificIntensity) const
    {
        EVector directv = directDissociationIntegralv(specificIntensity);
        EVector solomonv = spontaneousDissociationSinkv();
        return directv + solomonv;
    }

    EVector BigH2Model::directDissociationIntegralv(const Spectrum& specificIntensity, bool heatRate) const
    {
        // TODO: this can take up quite some time. If we assume that the specific intensity is
        // always constant during the lifetime of a BigH2Model, the we can just calculate this
        // once, at creation. In practice, one run of the gas module will then only require one
        // call of this function.

        EVector result{EVector::Zero(_h2Data->numLv())};

        // For each level that has cross section data
        for (size_t iLv : _h2Data->levelsWithCrossSectionv())
        {
            // For each cross section
            for (const Spectrum& cs : _h2Data->directDissociationCrossSections(iLv))
            {
                // We will integrate over (part of, in case the input spectrum is not
                // wide enough) the grid for the cross section
                const Array& cs_nuv = cs.frequencyv();

                // Usable integration range
                double minFreq = std::max(cs.freqMin(), specificIntensity.freqMin());
                double maxFreq = std::min(cs.freqMax(), specificIntensity.freqMax());

                // Integration lower bound: Index right of the minimum frequency
                size_t iNuMin = TemplatedUtils::index(minFreq, cs_nuv);
                // Index left of the minimum frequency
                if (iNuMin > 0) iNuMin--;

                // Integration upper bound: Index right of the maximum frequency
                size_t iNuMax = TemplatedUtils::index(maxFreq, cs_nuv);

                // Integrand: flux * sigma
                // Start with cross section [cm-2]
                Array sigmaFv{cs.valuev()};
                for (size_t iNu = iNuMin; iNu <= iNuMax; iNu++)
                {
                    // Multiply with photon flux density [s-1 cm-2 Hz-1]: F_nu = 4pi
                    // I_nu / h nu. (As always constant factors are applied after
                    // integrating.)
                    double nu = cs_nuv[iNu];
                    sigmaFv[iNu] *= specificIntensity.evaluate(nu) / nu;

                    // If we are calculating the heating rate (erg s-1) instead of
                    // the number rate (s-1), multiply with the energy minus the
                    // threshold (i.e. the lowest frequency of the grid for that
                    // specific cross section)
                    if (heatRate) sigmaFv[iNu] *= Constant::PLANCK * (nu - cs_nuv[0]);
                }

                // Integrate to total number of dissociations (s-1)
                result(iLv) += Constant::FPI / Constant::PLANCK
                               * TemplatedUtils::integrate<double>(cs_nuv, sigmaFv, iNuMin, iNuMax);
            }
        }
        return result;
    }

    const EVector& BigH2Model::spontaneousDissociationSinkv() const { return _h2Data->dissociationProbabilityv(); }
}
