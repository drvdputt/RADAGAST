#include "BigH2Model.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "H2Data.hpp"
#include "LevelSolver.hpp"
#include <iomanip>
#include <sstream>

namespace GasModule
{
    BigH2Model::BigH2Model(const H2Data* h2Data, const Spectrum* specificIntensity)
        : _h2Data{h2Data}, _specificIntensity{specificIntensity}, _levelSolution(_h2Data)
    {
        // since the radiation field stays constant, precalculate these integrals
        _directDissociationRatev = directDissociationIntegralv();
        _directDissociationHeatv = directDissociationIntegralv(true);
    }

    void BigH2Model::solve(double n, const CollisionParameters& cp, double h2form)
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

        _tvv = _h2Data->totalTransitionRatesvv(*_specificIntensity, cp, &_cvv, &_bpvv);
        _levelSolution.setCvv(_cvv);

        // TODO: Use better formation pumping recipe
        EVector sourcev = EVector::Zero(_h2Data->numLv());
        sourcev.head(_h2Data->startOfExcitedIndices()) = _h2Data->formationDistribution();
        sourcev *= h2form / sourcev.sum();

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

        // direct + solomon
        EVector sinkv = _directDissociationRatev + _h2Data->dissociationProbabilityv();
        // There are no transitions between electronically excited levels. Get this index here,
        // so that the solver called below can simplify the calculation using this assumption.
        int fullyConnectedCutoff = _h2Data->startOfExcitedIndices();
        EVector newNv = LevelSolver::statisticalEquilibrium_iterative(n, _tvv, sourcev, sinkv, _levelSolution.nv(),
                                                                      fullyConnectedCutoff);
        _levelSolution.setNv(newNv);
    }

    double BigH2Model::dissociationRate() const
    {
#ifdef STERNBERG2014
        // TODO: move this to simple h2 model
        // See 2014-Sternberg eq 3
        auto iv = *_specificIntensity.valuev();
        auto nuv = *_specificIntensity.frequencyv();

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
        // use fractional abundances here to get dissociation rate in [s-1]
        EVector fv = _levelSolution.fv();
        double directFractional = _directDissociationRatev.dot(fv);
        double solomonFractional = _h2Data->dissociationProbabilityv().dot(fv);
        DEBUG("Dissociation: direct rate:" << directFractional << " solomon rate: " << solomonFractional << '\n');
        return directFractional + solomonFractional;
#endif
    }

    double BigH2Model::dissociationHeating() const
    {
        // Fraction that dissociates per second * kinetic energy per dissociation * density of
        // level population = heating power density
        auto p = _h2Data->dissociationProbabilityv().array();
        auto k = _h2Data->dissociationKineticEnergyv().array();
        double solomonHeat = _levelSolution.nv().dot((p * k).matrix());
        double directHeat = _levelSolution.nv().dot(_directDissociationHeatv);
        DEBUG("Dissociation heat: direct heat: " << directHeat << " solomon heat:" << solomonHeat << '\n');
        return solomonHeat + directHeat;
    }

    double BigH2Model::netHeating() const
    {
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

        // Then add the dissociation cross section contribution by each level
        for (int i : _h2Data->levelsWithCrossSectionv())
        {
            for (const Spectrum& cs : _h2Data->directDissociationCrossSections(i))
                totalOpv += _levelSolution.nv()(i) * cs.binned(oFrequencyv);
        }
        return totalOpv;
    }

    void BigH2Model::extraDiagnostics(GasDiagnostics& gd) const
    {
        auto levelLabel = [](const H2Data::H2Level& level) {
            std::stringstream label;
            label << "E" << static_cast<int>(level.eState()) << std::setfill('0') << " J" << std::setw(2) << level.j()
                  << " v" << std::setw(2) << level.v();
            return label.str();
        };

        EVector fv = _levelSolution.fv();

        gd.setUserValue("H2 contdiss", _directDissociationRatev.dot(fv));
        for (int i : _h2Data->levelsWithCrossSectionv())
            gd.setUserValue("H2 contdiss " + levelLabel(_h2Data->level(i)), _directDissociationRatev[i] * fv[i]);

        const EVector& solomonDissv = _h2Data->dissociationProbabilityv();
        gd.setUserValue("H2 solomon", solomonDissv.dot(fv));
        for (int i = 0; i < _h2Data->numLv(); i++)
        {
            double dissContribution = solomonDissv[i] * fv[i];
            if (dissContribution > 0) gd.setUserValue("H2 solomon " + levelLabel(_h2Data->level(i)), dissContribution);
        }
    }

    EVector BigH2Model::directDissociationIntegralv(bool heatRate) const
    {
        EVector result{EVector::Zero(_h2Data->numLv())};

        // For each level that has cross section data
        for (int i : _h2Data->levelsWithCrossSectionv())
        {
            // For each cross section
            for (const Spectrum& cs : _h2Data->directDissociationCrossSections(i))
            {
                // We will integrate over (part of, in case the input spectrum is not
                // wide enough) the grid for the cross section
                const Array& cs_nuv = cs.frequencyv();

                // Usable integration range
                double minFreq = std::max(cs.freqMin(), _specificIntensity->freqMin());
                double maxFreq = std::min(cs.freqMax(), _specificIntensity->freqMax());

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
                    sigmaFv[iNu] *= _specificIntensity->evaluate(nu) / nu;

                    // If we are calculating the heating rate (erg s-1) instead of
                    // the number rate (s-1), multiply with the energy minus the
                    // threshold (i.e. the lowest frequency of the grid for that
                    // specific cross section)
                    if (heatRate) sigmaFv[iNu] *= Constant::PLANCK * (nu - cs_nuv[0]);
                }

                // Integrate to total number of dissociations (s-1)
                result(i) += Constant::FPI / Constant::PLANCK
                             * TemplatedUtils::integrate<double>(cs_nuv, sigmaFv, iNuMin, iNuMax);
            }
        }
        return result;
    }
}
