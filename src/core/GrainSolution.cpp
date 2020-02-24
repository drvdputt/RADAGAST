#include "GrainSolution.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GrainH2Formation.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "GrainPhotoelectricData.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "TemplatedUtils.hpp"

GrainSolution::GrainSolution(const GasModule::GrainPopulation* population)
    : _population{population},
      _chargeDistributionv(population->numSizes()), _newTemperaturev{population->initialTemperaturev()},
      _h2Heatv(population->numSizes())
{
    if (population->photoelectricData())
        _photoelectricCalculator = population->photoelectricData()->makeCalculator(population->sizev());
}

void GrainSolution::recalculate(const Spectrum* specificIntensity, double T, const SpeciesVector& sv)
{
    _specificIntensity = specificIntensity;
    _photoelectricLocals = GrainPhotoelectricCalculator::Locals(
        specificIntensity, T, sv.ne(), {-1, 1, 0, 0}, {sv.ne(), sv.np(), sv.nH(), sv.nH2()},
        {Constant::ELECTRONMASS, Constant::PROTONMASS, Constant::HMASS_CGS, 2 * Constant::HMASS_CGS});

    if (_population->h2formationData())
        _h2Heatv = _population->h2formationData()->surfaceH2FormationHeatPerSize(_population->sizev(), _newTemperaturev,
                                                                                 T, sv.nH());

    // I think the charge distributions are best calculated first, since they affect the the
    // temperatures, but not vice versa.
    recalculateChargeDistributions();
    recalculateTemperatures();
}

void GrainSolution::recalculateTemperatures()
{
    // Heat going into grains due to h2 formation on their surfaces
    Array extraGrainHeatv = _h2Heatv;

    if (_photoelectricCalculator && Options::cooling_gasGrainCollisions)
    {
        for (int i = 0; i < _population->numSizes(); i++)
        {
            // Heat going into grain due to collisions
            extraGrainHeatv[i] += _photoelectricCalculator->gasGrainCollisionCooling(
                i, _photoelectricLocals, _chargeDistributionv[i], _newTemperaturev[i], true);

            // Reduction of heating due to ejection of photoelectrons. TODO: cache this if slow (is
            // calculated in photoelectricGasHeating too)
            double grainPhotoPerSize = _photoelectricCalculator->heatingRateA(
                i, _photoelectricLocals, _population->qAbsv(i), _chargeDistributionv[i]);
            extraGrainHeatv[i] -= grainPhotoPerSize;
        }
    }

    DEBUG("New temps for grains:");
    for (size_t i = 0; i < _population->numSizes(); i++)
    {
        double cross = Constant::PI * _population->size(i) * _population->size(i);
        // work with the total absorption here (hence we multiply by 4pi and pi a^2). This makes it
        // more straightforward to factor in the 'extra' heating contributions.
        double absorption = Constant::FPI * cross
                            * TemplatedUtils::integrate<double, Array, Array>(
                                _specificIntensity->frequencyv(), _population->qAbsv(i) * _specificIntensity->valuev());

        auto heating = [&](double T) -> int {
            double bbEmission = 0;
            if (T > 0.) bbEmission = _population->totalThermalEmission(i, T);
            if (bbEmission < absorption + extraGrainHeatv[i]) return 1;
            if (bbEmission > absorption + extraGrainHeatv[i])
                return -1;
            else
                return 0;
        };
        _newTemperaturev[i] = TemplatedUtils::binaryIntervalSearch<double>(
            heating, 30., 1.e-3, Options::grainsolution_maxGrainTemp, Options::grainsolution_minGrainTemp);
        DEBUG(' ' << i << " " << _newTemperaturev[i]);
    }
    DEBUG('\n');
}

void GrainSolution::recalculateChargeDistributions()
{
    if (!_photoelectricCalculator) return;

    DEBUG("grain average charge:");
    for (int i = 0; i < _population->numSizes(); i++)
    {
        _photoelectricCalculator->calculateChargeDistribution(i, _photoelectricLocals, _population->qAbsv(i),
                                                              _chargeDistributionv[i]);
        DEBUG(' ' << i << ' ' << _chargeDistributionv[i].average());
    }
    DEBUG('\n');
}

double GrainSolution::photoelectricGasHeating()
{
    if (!_photoelectricCalculator) return 0;

    DEBUG("grain heat contributions: ");
    double total = 0;
    for (int i = 0; i < _population->numSizes(); i++)
    {
        double contribution = _population->density(i)
                 * _photoelectricCalculator->heatingRateA(i, _photoelectricLocals, _population->qAbsv(i),
                                                          _chargeDistributionv[i]);
        total += contribution;

        DEBUG(' ' << i << ' ' << contribution);

        // The net heating rate (eq 41 without denominator)
        if (Options::grainphotoelectriceffect_recombinationCooling)
            total -=
                _population->density(i)
                * _photoelectricCalculator->recombinationCoolingRate(i, _photoelectricLocals, _chargeDistributionv[i]);
    }
    DEBUG('\n');
    return total;
}

double GrainSolution::collisionalGasCooling() const
{
    if (!_photoelectricCalculator) return 0;

    double total = 0;
    for (int i = 0; i < _population->numSizes(); i++)
        total += _population->density(i)
                 * _photoelectricCalculator->gasGrainCollisionCooling(i, _photoelectricLocals, _chargeDistributionv[i],
                                                                      _newTemperaturev[i], false);
    return total;
}

double GrainSolution::surfaceH2FormationRateCoeff(double Tgas) const
{
    if (!_population->h2formationData()) return 0;

    return _population->h2formationData()->surfaceH2FormationRateCoeff(_population->sizev(), _newTemperaturev,
                                                                       _population->densityv(), Tgas);
}
