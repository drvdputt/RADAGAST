#ifndef CORE_GLOVERABEL08_HPP
#define CORE_GLOVERABEL08_HPP

#include <vector>

/** Coefficients and equations from Glover and Abel 2008, MNRAS, 388, 1627. Collisions with He are
    also available, but are not implemented. */
namespace GloverAbel08
{
    // Table 1. Fitting coefficients for the cooling rate of ortho-H2 excited by collisions with
    // atomic hydrogen.
    const std::vector<double> ortho_H_100to1000K_av = {-24.330855, 4.4404496, -4.0460989,
                                                       -1.1390725, 9.8094223, 8.6273872 /*, 0.0, 0.0*/};
    const std::vector<double> ortho_H_1000to6000K_av = {-24.329086, 4.6105087, -3.9505350, 12.363818,
                                                        -32.403165, 48.853562, -38.542008, 12.066770};

    // Table 2. Fitting coefficients for the cooling rate of para-H2 excited by collisions with
    // atomic hydrogen.
    const std::vector<double> para_H_100to1000K_av = {-24.216387, 3.3237480,  -11.642384,
                                                      -35.553366, -35.105689, -10.922078 /*, 0.0, 0.0*/};
    const std::vector<double> para_H_1000to6000K_av = {-24.216387, 4.2046488,   -1.3155285, -1.6552763,
                                                       4.1780102,  -0.56949697, -3.3824407, 1.0904027};

    // Table 3. Fitting coefficients for the cooling rate of para-H2 excited by collisions with H2.
    const std::vector<double> para_para_av = {-23.889798, 1.8550774, -0.55593388, 0.28429361, -0.20581113, 0.13112378};
    const std::vector<double> para_ortho_av = {-23.748534, 1.76676480,  -0.58634325,
                                               0.31074159, -0.17455629, 0.18530758};

    // Table 4. Fitting coefficients for the cooling rate of ortho-H2 excited by collisions with H2.
    const std::vector<double> ortho_para_av = {-24.126177, 2.3258217, -1.0082491, 0.54823768, -0.33679759, 0.20771406};
    const std::vector<double> ortho_ortho_av = {-24.020047, 2.2687566,   -1.0200304,
                                                0.83561432, -0.40772247, 0.096025713};

    // Table 6. Fitting coefficients for the cooling rate of H2 excited by collisions with protons.
    const std::vector<double> para_p_av = {-21.757160, 1.3998367, -0.37209530, 0.061554519, -0.37238286, 0.23314157};
    const std::vector<double> ortho_p_av = {-21.706641, 1.3901283, -0.34993699, 0.075402398, -0.23170723, 0.068938876};

    // Table 7. Fitting coefficients for the cooling rate of H2 excited by collisions with
    // electrons.
    const std::vector<double> para_e_below1000K_av = {-22.817869, 0.95653474, 0.79283462,
                                                      0.56811779, 0.27895033, 0.056049813};
    const std::vector<double> para_e_above1000K_av = {-22.817869, 0.66916141, 7.1191428,
                                                      -11.176835, 7.0467275,  -1.6471816};
    const std::vector<double> ortho_e_av = {-21.703215, 0.76059565, 0.50644890, 0.050371349, -0.10372467, -0.035709409};

    /** Cooling rate due to ortho H2 excited by H collisions. [erg cm3 s-1] */
    double coolOrthoH(double T);
    /** Cooling rate due to para H2 excited by H collisions. [erg cm3 s-1] */
    double coolParaH(double T);

    /** Cooling rate due to para H2 excited by para H2 collisions. [erg cm3 s-1] */
    double coolParaPara(double T);
    /** Cooling rate due to para H2 excited by ortho H2 collisions. [erg cm3 s-1] */
    double coolParaOrtho(double T);
    /** Cooling rate due to ortho H2 excited by para H2 collisions. [erg cm3 s-1] */
    double coolOrthoPara(double T);
    /** Cooling rate due to ortho H2 excited by ortho H2 collisions. [erg cm3 s-1] */
    double coolOrthoOrtho(double T);
    /** Cooling rate due to H2 (either ortho or para) excited by collisions with H2 (either ortho or
        para). The second argument should be one of the four H2 collision coefficient sets. */
    double coolH2H2Polynomial(double T, const std::vector<double>& coefficients);

    /** Cooling rate due to para H2 excited by p+ collisions. [erg cm3 s-1] */
    double coolParaProton(double T);
    /** Cooling rate due to ortho H2 excited by p+ collisions. [erg cm3 s-1] */
    double coolOrthoProton(double T);
    /** Cooling rate due to H2 (either ortho or para) excited by collisions with protons. The second
        argument should be one of the two p+ collision coefficient sets. */
    double coolProtonPolynomial(double T, const std::vector<double>& coefficients);

    /** Cooling rate due to ortho-para conversions induced by p+ collisions. Should only matter when
        out of thermodynamic equilibrium. [erg cm3 s-1] */
    double coolOrthoParaConversionProton(double T, double orthoFrac);

    /** Cooling rate due to para H2 excited by e- collisions. [erg cm3 s-1] */
    double coolParaElectron(double T);
    /** Cooling rate due to ortho H2 excited by e- collisions. [erg cm3 s-1] */
    double coolOrthoElectron(double T);

    /** The x_k argument is the x parameter of equation 36 divided by the boltzman constant */
    double coolElectronPolynomial(double T, double x_k, const std::vector<double>& coefficients);

}  // namespace GloverAbel08

#endif  // CORE_GLOVERABEL08_HPP
