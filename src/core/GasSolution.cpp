#include "GasSolution.hpp"
#include "GasInterfaceImpl.hpp"
#include "GrainInterface.hpp"
#include "GrainType.hpp"
#include "IonizationBalance.hpp"
#include "Options.hpp"

Array GasSolution::emisivityv(const Array& eFrequencyv) const
{
	Array lineEmv(eFrequencyv.size());
	lineEmv += _hSolution->emissivityv(eFrequencyv);
	lineEmv += _h2Solution->emissivityv(eFrequencyv);

	Array contEmCoeffv(eFrequencyv.size());
	contEmCoeffv += _gasInterfaceImpl->radiativeRecombinationEmissivityv(_t, eFrequencyv);
	contEmCoeffv += _gasInterfaceImpl->freeFreeEmissivityv(_t, eFrequencyv);

	return lineEmv + (np() * ne() / Constant::FPI) * contEmCoeffv;
}

Array GasSolution::opacityv(const Array& oFrequencyv) const
{
	size_t numFreq = oFrequencyv.size();
	Array lineOpv(numFreq);
	lineOpv += _hSolution->opacityv(oFrequencyv);
	lineOpv += _h2Solution->opacityv(oFrequencyv);

	Array contOpv = np() * ne() * _gasInterfaceImpl->freeFreeOpacityv(_t, oFrequencyv);

	Array totalOpv(numFreq);
	for (size_t i = 0; i < numFreq; i++)
	{
		// TODO: this should actually be the average over the cross section for
		// this frequency bin
		double ionizOp_iFreq = nH() * Ionization::crossSection(oFrequencyv[i]);
		totalOpv[i] = ionizOp_iFreq + contOpv[i] + lineOpv[i];
	}
	return totalOpv;
}

double GasSolution::cooling() const
{
	double lineCool = -(_hSolution->netHeating() + _h2Solution->netHeating());
	double freefreeCool = _gasInterfaceImpl->freeFreeCool(np() * ne(), _t);
	double hRecCool = Ionization::cooling(nH(), np(), ne(), _t);
	return lineCool + freefreeCool + hRecCool;
}

double GasSolution::heating() const
{
	// double freefreeHeat = _freeFree->heating(np(s) * ne(s), s.T, s.specificIntensity);
	// TODO: decide whether to keep the above, as it is negligible in any case I can imagine
	double hPhotoIonHeat = Ionization::heating(np(), ne(), _t, _specificIntensity);
	double dissHeat = _h2Solution->dissociationHeating(_specificIntensity);
	double grainHeat = grainHeating();
	return hPhotoIonHeat + dissHeat + grainHeat;
}

double GasSolution::grainHeating(double* photoHeat, double* collCool) const
{
	double grainPhotoelectricHeating = 0;
	double gasGrainCooling = 0;

	// Specify the environment parameters
	GrainPhotoelectricEffect::Environment env(_specificIntensity, _t, ne(), np(),
	                                          {-1, 1, 0, 0}, {ne(), np(), nH(), nH2()},
	                                          {Constant::ELECTRONMASS, Constant::PROTONMASS,
	                                           Constant::HMASS_CGS,
	                                           2 * Constant::HMASS_CGS});

	size_t numPop = _grainInterface.numPopulations();
	for (size_t i = 0; i < numPop; i++)
	{
		const GasModule::GrainInterface::Population* pop =
		                _grainInterface.population(i);
		const GrainType* type = pop->type();
		if (type->heatingAvailable())
		{
			/* Choose the correct parameters for the photoelectric effect based on
			   the type (a.k.a. composition) of the Population. */
			GrainPhotoelectricEffect gpe(*type);

			size_t numSizes = pop->numSizes();
			for (int m = 0; m < numSizes; m++)
			{
				double a = pop->size(m);
				const Array& qAbsv = pop->qAbsv(m);
				double nd = pop->density(m);
				auto cd = gpe.calculateChargeDistribution(a, env, qAbsv);
				grainPhotoelectricHeating +=
				                nd * gpe.heatingRateA(a, env, qAbsv, cd);
				if (Options::cooling_gasGrainCollisions)
					gasGrainCooling += nd *
					                   gpe.gasGrainCollisionCooling(
					                                   a, env, cd,
					                                   pop->temperature(m),
					                                   false);
			}
		}
	}
	if (photoHeat != nullptr)
		*photoHeat = grainPhotoelectricHeating;
	if (collCool != nullptr)
		*collCool = gasGrainCooling;
	return grainPhotoelectricHeating - gasGrainCooling;
}
