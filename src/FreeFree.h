#ifndef _SRC_FREEFREE_H_
#define _SRC_FREEFREE_H_

#include "Table.h"

#include <string>
#include <vector>

class FreeFree
{
public:
	FreeFree(const Array& frequencyv);

private:
	void readData(const std::string& filename);

public:
	/* Calculate the mission coefficient for the free-free continuum for all frequencies. The units are
	[density^-1][power]/[frequency interval] cm^3 erg / s / cm. The emissivity
	([power][density]/[frequency interval]) can be obtained by multiplying this value with ne*np / 4pi. The
	contributions are added to the current contents of gamma_nu. */
	void addEmissionCoefficientv(double T, Array& gamma_nuv) const;

	/* Calculate the opacity coefficient for the free-free continuum for all frequencies. The units are
	 [density^-2][length^-1]. Multiplying with ne*np will give the opacity in [length-1]. The contirbutiona are added
	 to the current contents of */
	void addOpacityCoefficientv(double T, Array& opCoeffv) const;

	double gauntFactor(double logu, double logg2) const;

private:
	const Array& _frequencyv;

	/* The ((natural!) log of the) data read from the file, indexed on (u, gamma^2) */
	Table<2> _fileGauntFactorvv;
	/* The ((base 10!) log of the) minima and maxima defining the u - gamma^2 grid, as well as the step
	 size in log space. */
	double _loggamma2Min, _loguMin, _logStep, _loggamma2Max, _loguMax;
};

#endif /* _SRC_FREEFREE_H_ */
