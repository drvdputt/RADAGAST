#include "GrainSolution.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "GrainPhotoelectricData.hpp"
#include "GrainPhotoelectricEffect.hpp"
#include "GrainPopulation.hpp"
#include "SpecialFunctions.hpp"
#include "Spectrum.hpp"
#include "TemplatedUtils.hpp"

GrainSolution::GrainSolution(const GrainPopulation& population)
                : _population{&population}, _chargeDistributionv(population.numSizes()),
                  _newTemperaturev{population.temperaturev()}
{
	if (population.photoelectricData())
		_photoelectricCalculator = population.photoelectricData()->makeCalculator(
		                population.sizev());
}

void GrainSolution::recalculateTemperatures(const Spectrum& specificIntensity,
                                            Array otherGrainHeat)
{
	const Array& frequencyv = specificIntensity.frequencyv();
	const auto& qAbsvv = _population->qAbsvv();

	for (size_t i = 0; i < _population->numSizes(); i++)
	{
		double cross = Constant::PI * _population->size(i) * _population->size(i);
		double absorption =
		                cross * TemplatedUtils::integrate<double, Array, Array>(
		                                        frequencyv,
		                                        qAbsvv[i] * specificIntensity.valuev());

		auto heating = [&](double T) -> int {
			Array blackbodyIntegrandv(frequencyv.size());
			for (size_t j = 0; j < frequencyv.size(); j++)
				blackbodyIntegrandv[j] =
				                qAbsvv[i][j] *
				                SpecialFunctions::planck(frequencyv[j], T);

			double bbEmission = 0;
			if (T > 0.)
				bbEmission = cross *
				             TemplatedUtils::integrate<double, Array, Array>(
				                             frequencyv, blackbodyIntegrandv);
			if (bbEmission < absorption + otherGrainHeat[i])
				return 1;
			if (bbEmission > absorption + otherGrainHeat[i])
				return -1;
			else
				return 0;
		};
		_newTemperaturev[i] = TemplatedUtils::binaryIntervalSearch<double>(
		                heating, 30., 1.e-3, 300, 1.);
		DEBUG("New temp for grain " << i << " " << _newTemperaturev[i] << " K\n");
	}
}

void GrainSolution::recalculateChargeDistributions(
                const GrainPhotoelectricCalculator::Environment& env)
{
	if (!_photoelectricCalculator)
		return;

	for (int i = 0; i < _population->numSizes(); i++)
		_chargeDistributionv[i] = _photoelectricCalculator->calculateChargeDistribution(
		                i, env, _population->qAbsv(i));
}

double GrainSolution::photoelectricGasHeating(
                const GrainPhotoelectricCalculator::Environment& env) const
{
	if (!_photoelectricCalculator)
		return 0;

	double total = 0;

	for (int i = 0; i < _population->numSizes(); i++)
	{
		const Array& qAbsv = _population->qAbsv(i);
		double nd = _population->density(i);
		total += nd * _photoelectricCalculator->heatingRateA(i, env, qAbsv,
		                                                     _chargeDistributionv[i]);
	}
	return total;
}

double
GrainSolution::collisionalGasCooling(const GrainPhotoelectricCalculator::Environment& env) const
{
	if (!_photoelectricCalculator)
		return 0;

	double total = 0;
	for (int i = 0; i < _population->numSizes(); i++)
	{
		double nd = _population->density(i);
		total += nd * _photoelectricCalculator->gasGrainCollisionCooling(
		                              i, env, _chargeDistributionv[i],
		                              _population->temperature(i), false);
	}
	return total;
}
