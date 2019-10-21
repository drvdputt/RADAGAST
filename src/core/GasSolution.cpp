#include "GasSolution.hpp"
#include "GasInterfaceImpl.hpp"
#include "IonizationBalance.hpp"

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
	lineOp += _hSolution->opacityv(oFrequencyv);
	lineOp += _h2Solution->opacityv(oFrequencyv);

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

