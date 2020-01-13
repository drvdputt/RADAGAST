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

GrainSolution::GrainSolution(const GrainPopulation& population)
    : _population{&population}, _chargeDistributionv(population.numSizes()), _newTemperaturev{population.temperaturev()}
{
    if (population.photoelectricData())
        _photoelectricCalculator = population.photoelectricData()->makeCalculator(population.sizev());
}

void GrainSolution::recalculateTemperatures(const GrainPhotoelectricCalculator::Environment& env)
{
    Array extraGrainHeatPerSize(_population->numSizes());

    if (_photoelectricCalculator && Options::cooling_gasGrainCollisions)
    {
        for (int i = 0; i < _population->numSizes(); i++)
        {
            extraGrainHeatPerSize[i] += _photoelectricCalculator->gasGrainCollisionCooling(
                _population->size(i), env, _chargeDistributionv[i], _newTemperaturev[i], true);
            // cout << "extra grain heat " << m << " " << grainHeatPerSizev[m] << '\n';

            // grainPhotoPerSizev[m] = gpe.heatingRateA(pop.size(m), env, pop.qAbsv(m), cd);
            // cout << "- photo heat " << m << " " << grainPhotoPerSizev[m] << '\n';
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
        _chargeDistributionv[i] = _photoelectricCalculator->calculateChargeDistribution(i, env, _population->qAbsv(i));
}

double GrainSolution::photoelectricGasHeating(const GrainPhotoelectricCalculator::Environment& env) const
{
    if (!_photoelectricCalculator) return 0;

    double total = 0;

    for (int i = 0; i < _population->numSizes(); i++)
    {
        const Array& qAbsv = _population->qAbsv(i);
        double nd = _population->density(i);
        total += nd * _photoelectricCalculator->heatingRateA(i, env, qAbsv, _chargeDistributionv[i]);
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
