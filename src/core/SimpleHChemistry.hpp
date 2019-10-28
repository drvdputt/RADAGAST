#ifndef CORE_SIMPLEHCHEMISTRY_HPP
#define CORE_SIMPLEHCHEMISTRY_HPP

#include "Chemistry.hpp"

class Spectrum;

class SimpleHChemistry : public Chemistry
{
public:
	SimpleHChemistry();
	EVector rateCoeffv(double T, const Spectrum& specificIntensity,
	                   double kDissFromH2Levels, double kH2FormationGrain) const;
};

#endif // CORE_SIMPLEHCHEMISTRY_HPP
