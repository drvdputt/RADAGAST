#ifndef _FREEFREE_H_
#define _FREEFREE_H_

#include "Table.h"

#include <string>
#include <vector>

class FreeFree
{
public:
	FreeFree(const Array& frequencyv);

private:
	void readFullData(const std::string& file);
	void readIntegratedData(const std::string& file);
	double gauntFactor(double logu, double logg2) const;
	double integratedGauntFactor(double logg2) const;

public:
	/* Calculate the emission coefficient for the free-free continuum for all frequencies. The
	units are [density^-1][power]/[frequency interval] cm^3 erg / s / cm. The emissivity
	([power][density]/[frequency interval]) can be obtained by multiplying this value with ne*np
	/ 4pi. The contributions are added to the current contents of gamma_nu. */
	void addEmissionCoefficientv(double T, Array& gamma_nuv) const;

	/* Calculate the opacity coefficient for the free-free continuum for all frequencies. The
	 units are [density^-2][length^-1]. Multiplying with ne*np will give the opacity in
	 [length-1]. The contributions are added to the current contents of */
	void addOpacityCoefficientv(double T, Array& opCoeffv) const;

	/* Calculate the heating due to direct absorption of photons the free electrons moving in
	 * the field of the free protons. */
	double heating(double np_ne, double T, const Array& specificIntensityv) const;

	/* Calculate the cooling due to Bremsstrahlung, given the product of the electron and proton
	 * densities, and their kinetic temperature. */
	double cooling(double np_ne, double T) const;

private:
	Array _frequencyv;

	/* The ((natural!) log of the) data read from the file, indexed on (u, gamma^2) */
	Table<2> _fileGauntFactorvv;
	/* The ((base 10!) log of the) minima and maxima defining the u - gamma^2 grid, as well as
	 the step size in log space. */
	double _loggamma2Min{0}, _loguMin{0}, _logStep{0}, _loggamma2Max{0}, _loguMax{0};

	/* The natural log of the data read from the file containing the integrated Gaunt factor,
	 * indexed on (gamma^2). Note that there is a different number of data points here in the
	 * gamma^2 direction. I'm not really sure why we store the gamma^2 values here. Right now
	 * I'm just going to adapt the implementation from interpolate4.c */
	Array _loggamma2_integrated, _fileGauntFactorv_integrated;

	/* The values specifying the log10(gamma^2) grid */
	double _loggamma2Min_integrated{0}, _logStep_integrated{0}, _loggamma2Max_integrated{0};
};

#endif /* _FREEFREE_H_ */
