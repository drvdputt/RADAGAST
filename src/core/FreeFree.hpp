#ifndef CORE_FREEFREE_HPP
#define CORE_FREEFREE_HPP

#include "Spectrum.hpp"
#include "Table.hpp"
#include <string>

namespace RADAGAST
{
    /** This class provides implementations for quantities related to the free-free processes in a
        hydrogen plasma. It is mainly based on the data for the free-free Gaunt factor that I got
        from a 2014 paper from van Hoof et al. (MNRAS 444 420). */
    class FreeFree
    {
    public:
        /** All the data is read in during construction. The given frequency grid will be used for
            the output. */
        FreeFree();

    private:
        /** Reads in gauntff_merged_Z01.dat. */
        void readFullData();

        /** Reads in gauntff_freqint_Z01.dat. */
        void readIntegratedData();

        /** Returns the gaunt factor, interpolated for a specific log(u) and log(gamma^2). */
        double gauntFactor(double logu, double logg2) const;

        /** Return the frequency-integrated guant factor, interpolated for a specific
            log(gamma^2). */
        double integratedGauntFactor(double logg2) const;

    public:
        /** Calculate the emission coefficient for the free-free continuum for all frequencies. The
            units are [density^-1][power]/[frequency interval] cm^3 erg / s / cm. The emissivity
            ([power][density]/[frequency interval]) can be obtained by multiplying this value with
            ne*np / 4pi. The contributions are added to the current contents of gamma_nu. */
        void addEmissionCoefficientv(double T, const Array& eFrequencyv, Array& gamma_nuv) const;

        /** Calculate the opacity coefficient for the free-free continuum for all frequencies. The
            units are [density^-2][length^-1]. Multiplying with ne*np will give the opacity in
            [length-1]. The contributions are added to the current contents of opCoeffv in
            [cm-7]. */
        void addOpacityCoefficientv(double T, const Array& oFrequencyv, Array& opCoeffv) const;

        /** Calculate the opacity coefficient for the free-free continuum for a specific frequency
            and temperature. The units of the returned coefficient is [cm-7]: opacity [cm-1] per
            density squared [cm-6]. Multiplying with ne*np will give the opacity in [cm-1]. */
        double opacityCoefficient(double nu, double T) const;

        /** Calculate the cooling due to Bremsstrahlung, given the product of the electron and
            proton densities, and their kinetic temperature. [erg cm-3 s-1]. */
        double cooling(double np_ne, double T) const;

    private:
        /* The ((natural!) log of the) data read from the file, indexed on (u, gamma^2) */
        Table<2> _fileGauntFactorvv;
        /* The ((base 10!) log of the) minima and maxima defining the u - gamma^2 grid, as well as
           the step size in log space. */
        double _loggamma2Min{0}, _loguMin{0}, _logStep{0}, _loggamma2Max{0}, _loguMax{0};

        /* The natural log of the data read from the file containing the integrated Gaunt factor,
           indexed on (gamma^2). Note that there is a different number of data points here in the
           gamma^2 direction. I'm not really sure why we store the gamma^2 values here. Right now
           I'm just going to adapt the implementation from interpolate4.c */
        Array _loggamma2_integrated, _fileGauntFactorv_integrated;

        /* The values specifying the log10(gamma^2) grid */
        double _loggamma2Min_integrated{0}, _logStep_integrated{0}, _loggamma2Max_integrated{0};
    };
}
#endif  // CORE_FREEFREE_HPP
