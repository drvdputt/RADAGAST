#include "RadiationFieldTools.hpp"
#include "Constants.hpp"
#include "SpecialFunctions.hpp"

#include <vector>

Array RadiationFieldTools::generateSpecificIntensityv(const Array& frequencyv, double Tc,
                                                      double G0)
{
	Array I_nu(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		I_nu[iFreq] = SpecialFunctions::planck(frequencyv[iFreq], Tc);

	// Cut out the UV part
	size_t i = 0;
	size_t startUV, endUV;
	while (frequencyv[i] < Constant::LIGHT / (2400 * Constant::ANG_CM) &&
	       i < frequencyv.size())
		i++;
	startUV = i > 0 ? i - 1 : 0;
	while (frequencyv[i] < Constant::LIGHT / (912 * Constant::ANG_CM) &&
	       i < frequencyv.size())
		i++;
	endUV = i + 1;
	std::vector<double> frequenciesUV(begin(frequencyv) + startUV,
	                                  begin(frequencyv) + endUV);
	std::vector<double> isrfUV(begin(I_nu) + startUV, begin(I_nu) + endUV);

	// Integrate over the UV only
	double UVdensity = Constant::FPI / Constant::LIGHT *
	                   TemplatedUtils::integrate<double>(frequenciesUV, isrfUV);
	double currentG0 = UVdensity / Constant::HABING_DENS;

	// Rescale to _G0
	I_nu *= G0 / currentG0;

	std::vector<double> isrfUVbis(begin(I_nu) + startUV, begin(I_nu) + endUV);
	return I_nu;
}

Array RadiationFieldTools::freqToWavGrid(const Array& frequencyv)
{
	size_t numWav = frequencyv.size();
	Array wavelengthv(numWav);
	for (size_t iWav = 0; iWav < numWav; iWav++)
		wavelengthv[iWav] = Constant::LIGHT / frequencyv[numWav - 1 - iWav];
	return wavelengthv;
}

Array RadiationFieldTools::freqToWavSpecificIntensity(const Array& frequencyv,
                                                      const Array& specificIntensity_nu)
{
	Array I_lambda(frequencyv.size());
	for (size_t iFreq = 0; iFreq < frequencyv.size(); iFreq++)
		I_lambda[I_lambda.size() - iFreq - 1] = specificIntensity_nu[iFreq] *
		                                        frequencyv[iFreq] * frequencyv[iFreq] /
		                                        Constant::LIGHT;
	return I_lambda;
}
