#ifndef _FREEBOUND_H_
#define _FREEBOUND_H_

#include "Table.h"

#include <vector>

class FreeBound
{
public:
	/* Creates an object and reads in the free-bound continuum emission data. The data is
	   interpolated for every frequency on the grid, and stored for the temperatures contained
	   in the file. */
	FreeBound(const Array& frequencyv);

private:
	/* Function that reads the data. To be called in the constructor*/
	void readData(std::string file, std::vector<double>& fileFrequencyv,
	              std::vector<double>& fileThresholdv,
	              std::vector<double>& fileTemperaturev,
	              std::vector<std::vector<double>>& fileGammaDaggervv) const;

public:
	Array thresholdv() const { return _thresholdv; }

	/* Calculate the emission coefficient for the optical recombination continuum for
	   allfrequencies. The data is intepolated ad-hoc in the temperature direction; in the
	   frequency direction this data was already interpolated in the constructor. Returned in
	   units [density^-1][power]/[frequency interval] cm^3 erg / s / Hz. The emissivity
	   ([power][density]/[frequency interval]) can be obtained by multiplying this value with
	   ne_np / 4pi. The contribution at each frequency is added to the current contents of
	   gamma_nu */
	void addEmissionCoefficientv(double T, Array& gamma_nuv) const;

private:
	/* The frequency grid onto which the data will be interpolated */
	const Array& _frequencyv;

	/* Data to be loaded in constructor body. First index is for frequency, second for
	   temperature */
	Table<2> _gammaDaggervv;
	/* Vector containing the threshold frequencies, i.e. those of the lines starting with 1. Is
	   needed for applying equation 1 of Ercolano and Storey 2006 (MNRAS 372, 1875) */
	Array _thresholdv;

	/* The accompanying log-temperature grid */
	std::vector<double> _logTemperaturev;
};

#endif /* _FREEBOUND_H_ */
