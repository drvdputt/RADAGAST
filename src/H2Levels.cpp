#include "H2Levels.h"
#include "Constants.h"
#include "H2FromFiles.h"
#include "TemplatedUtils.h"

H2Levels::H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv)
                : NLevel(hff, frequencyv), _hff(hff)
{
}

H2Levels::~H2Levels() = default;

double H2Levels::dissociationRate(const NLevel::Solution& s, const Array& specificIntensityv) const
{
	// See 2014-Sternberg eq 3
	// F0 = integral 912 to 1108 Angstrom of Fnu(= 4pi Inu) with Inu in cm-2 s-1 Hz sr-1
	Array photonFluxv = Constant::FPI * specificIntensityv / frequencyv() / Constant::PLANCK;
	constexpr double freqLWmin{Constant::LIGHT / 1108 / Constant::ANG_CM};
	constexpr double freqLWmax{Constant::LIGHT / 912 / Constant::ANG_CM};
	size_t iLWmin{TemplatedUtils::index(freqLWmin, frequencyv())};
	size_t iLWmax{TemplatedUtils::index(freqLWmax, frequencyv())};
	double F0 = TemplatedUtils::integrate<double>(frequencyv(), photonFluxv, iLWmin, iLWmax);

	// eq 4 and 5
	double Iuv{F0 / 2.07e7};
	return 5.8e-11 * Iuv;
}

double H2Levels::dissociationHeating(const NLevel::Solution& s) const
{
	// TODO
	return 0;
}
