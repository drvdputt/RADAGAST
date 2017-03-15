#ifndef _FREEBOUND_H_
#define _FREEBOUND_H_

#include "ReadData.h"

class FreeBound
{
public:
	/* Creates an object and reads in the free-bound continuum emission data. The data is interpolated for
	/every frequency on the grid, and stored for the temperatures contained in the file. */
	FreeBound(const std::vector<double>& frequencyv);

	/* Calculate the emission coefficient for the optical recombination continuum for all frequencies.
	 The data is intepolated ad-hoc in the temperature direction; in the frequency direction this data was
	 already interpolated in the constructor. Returned in units [density^-1][power]/[frequency interval]
	 cm^3 erg / s / cm. The emissivity ([power][density]/[frequency interval]) can be obtained by
	 multiplying this value with ne_np */
	std::vector<double> emissionCoefficientv(double T) const;

private:
	// The frequency grid onto which the data will be interpolated
	const std::vector<double>& _frequencyv;

	/* Data to be loaded in constructor body. First index is for frequency, second for temperature */
	std::vector<std::vector<double>> _gammaDaggervv;
	/* Vector containing the threshold frequencies, i.e. those of the lines starting with 1. Is needed for
	/applying equation 1 of Ercolano and Storey 2006 (MNRAS 372, 1875) */
	std::vector<double> _thresholdv;

	// The accompanying log-temperature grid
	std::vector<double> _logTemperaturev;
};

#endif /* _FREEBOUND_H_ */
