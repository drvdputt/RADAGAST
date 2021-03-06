#ifndef CORE_FREEBOUND_HPP
#define CORE_FREEBOUND_HPP

#include "Table.hpp"
#include <string>
#include <vector>

namespace RADAGAST
{
    class FreeBound
    {
    public:
        /** Creates an object and reads in the free-bound continuum emission data. The data is
            ixnterpolated for every frequency on the grid, and stored for the temperatures contained
            in the file. */
        FreeBound();

    private:
        /** Function that reads the data. To be called in the constructor*/
        void readData(std::string file, std::vector<double>& fileFrequencyv, std::vector<double>& fileThresholdv,
                      std::vector<double>& fileTemperaturev, std::vector<std::vector<double>>& fileGammaDaggervv) const;

    public:
        const Array& thresholdv() const { return _thresholdv; }

        /** Calculate the emission coefficient for the optical recombination continuum for all
            frequencies. The data is intepolated ad-hoc in the temperature direction; in the
            frequency direction this data was already interpolated in the constructor. Returned
            in units [density^-1][power]/[frequency interval] cm^3 erg / s / Hz. The emissivity
            ([power][density]/[frequency interval][solid angle]) can be obtained by multiplying
            this value with ne_np / 4pi. The contribution at each frequency is added to the
            current contents of gamma_nu */
        void addEmissionCoefficientv(double T, const Array& eFrequencyv, Array& gamma_nuv) const;

    private:
        /** Emission coefficient for recombination to the ground state, using Maxwell
            distribution, Milne relation and ionization cross section. [erg s-1 Hz-1 cm3] */
        static double ionizingContinuumCoefficient(double T, double frequency);

        // Data to be loaded in constructor body. First index is for frequency, second for
        // temperature
        Array _frequencyv;
        Table<2> _gammaDaggervv;
        // Vector containing the threshold frequencies, i.e. those of the lines starting with 1. Is
        // needed for applying equation 1 of Ercolano and Storey 2006 (MNRAS 372, 1875)
        Array _thresholdv;

        // The accompanying log-temperature grid
        Array _logTemperaturev;
    };
}
#endif  // CORE_FREEBOUND_HPP
