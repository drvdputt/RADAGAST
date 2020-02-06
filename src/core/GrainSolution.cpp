#include "GrainSolution.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GrainH2Formation.hpp"
#include "GrainPhotoelectricCalculator.hpp"
#include "GrainPhotoelectricData.hpp"
#include "GrainPopulation.hpp"
#include "Options.hpp"
#include "SpecialFunctions.hpp"
#include "Spectrum.hpp"
#include "TemplatedUtils.hpp"

GrainSolution::GrainSolution(const GrainPopulation* population)
    : _population{population},
      _chargeDistributionv(population->sizev().size()), _newTemperaturev{population->temperaturev()}
{
    if (population->photoelectricData())
        _photoelectricCalculator = population->photoelectricData()->makeCalculator(population->sizev());
}

void GrainSolution::recalculateTemperatures(const GrainPhotoelectricCalculator::Environment& env)
{
    Array extraGrainHeatPerSize(_population->numSizes());

    if (_photoelectricCalculator && Options::cooling_gasGrainCollisions)
    {
        for (int i = 0; i < _population->numSizes(); i++)
        {
            extraGrainHeatPerSize[i] += _photoelectricCalculator->gasGrainCollisionCooling(
                i, env, _chargeDistributionv[i], _newTemperaturev[i], true);
            // cout << "extra grain heat " << m << " " << grainHeatPerSizev[m] << '\n';

            // TODO: cache this if slow (is calculated in photoelectricGasHeating too)
            double grainPhotoPerSize =
                _photoelectricCalculator->heatingRateA(i, env, _population->qAbsv(i), _chargeDistributionv[i]);
            extraGrainHeatPerSize[i] -= grainPhotoPerSize;
            // cout << "- photo heat " << m << " " << grainPhotoPerSize << '\n';
        }
    }

    if (_population->h2formationData())
    {
        // TODO use SpeciesVector here instead of hardcoded index
        extraGrainHeatPerSize += _population->h2formationData()->surfaceH2FormationHeatPerSize(
            _population->sizev(), _newTemperaturev, env._T, env._densityv[2]);
    }

    const Array& frequencyv = env._specificIntensity.frequencyv();
    const auto& qAbsvv = _population->qAbsvv();

    for (size_t i = 0; i < _population->numSizes(); i++)
    {
        double cross = Constant::PI * _population->size(i) * _population->size(i);
        double absorption =
            cross
            * TemplatedUtils::integrate<double, Array, Array>(frequencyv, qAbsvv[i] * env._specificIntensity.valuev());

        auto heating = [&](double T) -> int {
            double bbEmission = 0;
            if (T > 0.) bbEmission = _population->totalThermalEmission(i, T);
            if (bbEmission < absorption + extraGrainHeatPerSize[i]) return 1;
            if (bbEmission > absorption + extraGrainHeatPerSize[i])
                return -1;
            else
                return 0;
        };
        _newTemperaturev[i] = TemplatedUtils::binaryIntervalSearch<double>(
            heating, 30., 1.e-3, Options::grainsolution_maxGrainTemp, Options::grainsolution_minGrainTemp);
        DEBUG("New temp for grain " << i << " " << _newTemperaturev[i] << " K\n");
    }
}

void GrainSolution::recalculateChargeDistributions(const GrainPhotoelectricCalculator::Environment& env)
{
    if (!_photoelectricCalculator) return;

    for (int i = 0; i < _population->numSizes(); i++)
        _photoelectricCalculator->calculateChargeDistribution(i, env, _population->qAbsv(i), _chargeDistributionv[i]);
}

double GrainSolution::photoelectricGasHeating(const GrainPhotoelectricCalculator::Environment& env) const
{
    if (!_photoelectricCalculator) return 0;

    double total = 0;
    for (int i = 0; i < _population->numSizes(); i++)
    {
        total += _population->density(i)
                 * _photoelectricCalculator->heatingRateA(i, env, _population->qAbsv(i), _chargeDistributionv[i]);

        // The net heating rate (eq 41 without denominator)
        if (Options::grainphotoelectriceffect_recombinationCooling)
            total -= _population->density(i)
                     * _photoelectricCalculator->recombinationCoolingRate(i, env, _chargeDistributionv[i]);
    }
    return total;
}

double GrainSolution::collisionalGasCooling(const GrainPhotoelectricCalculator::Environment& env) const
{
    if (!_photoelectricCalculator) return 0;

    double total = 0;
    for (int i = 0; i < _population->numSizes(); i++)
    {
        double nd = _population->density(i);
        total += nd
                 * _photoelectricCalculator->gasGrainCollisionCooling(i, env, _chargeDistributionv[i],
                                                                      _population->temperature(i), false);
    }
    return total;
}

double GrainSolution::surfaceH2FormationRateCoeff(double Tgas) const
{
    if (!_population->h2formationData()) return 0;

    return _population->h2formationData()->surfaceH2FormationRateCoeff(_population->sizev(), _newTemperaturev,
                                                                       _population->densityv(), Tgas);
}
