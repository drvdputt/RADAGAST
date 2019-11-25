#ifndef CORE_GLOVERABEL08_HPP
#define CORE_GLOVERABEL08_HPP

#include <vector>

/** Coefficients and equations from Glover and Abel 2008, MNRAS, 388, 1627 */
namespace GloverAbel08
{
// Table 1. Fitting coefficients for the cooling rate of ortho-H2 excited by collisions with
// atomic hydrogen.
const std::vector<double> ortho_H_100to1000K_av = {-24.330855, 4.4404496,
                                                   -4.0460989, -1.1390725,
                                                   9.8094223,  8.6273872 /*, 0.0, 0.0*/};
const std::vector<double> ortho_H_1000to6000K_av = {-24.329086, 4.6105087,  -3.9505350,
                                                    12.363818,  -32.403165, 48.853562,
                                                    -38.542008, 12.066770};

// Table 2. Fitting coefficients for the cooling rate of para-H2 excited by collisions with
// atomic hydrogen.
const std::vector<double> para_H_100to1000K_av = {-24.216387, 3.3237480, -11.642384, -35.553366,
                                                  -35.105689, -10.922078 /*, 0.0, 0.0*/};
const std::vector<double> para_H_1000to6000K_av = {-24.216387, 4.2046488, -1.3155285,
                                                   -1.6552763, 4.1780102, -0.56949697,
                                                   -3.3824407, 1.0904027};

// Table 3. Fitting coefficients for the cooling rate of para-H2 excited by collisions with H2.
const std::vector<double> para_para_av = {-23.889798, 1.8550774,   -0.55593388,
                                          0.28429361, -0.20581113, 0.13112378};
const std::vector<double> para_ortho_av = {-23.748534, 1.76676480,  -0.58634325,
                                           0.31074159, -0.17455629, 0.18530758};

// Table 4. Fitting coefficients for the cooling rate of ortho-H2 excited by collisions with H2.
const std::vector<double> ortho_para_av = {-24.126177, 2.3258217,   -1.0082491,
                                           0.54823768, -0.33679759, 0.20771406};
const std::vector<double> ortho_ortho_av = {-24.020047, 2.2687566,   -1.0200304,
                                            0.83561432, -0.40772247, 0.096025713};

/** Cooling rate due to ortho H2 excited by H collisions. */
double coolOrthoH(double T);

/** Cooling rate due to para H2 excited by H collisions. */
double coolParaH(double T);

/** Cooling rate due to para H2 excited by para H2 collisions. */
double coolParaPara(double T);
/** Cooling rate due to para H2 excited by ortho H2 collisions. */
double coolParaOrtho(double T);
/** Cooling rate due to ortho H2 excited by para H2 collisions. */
double coolOrthoPara(double T);
/** Cooling rate due to ortho H2 excited by ortho H2 collisions. */
double coolOrthoOrtho(double T);

/** Cooling rate due to H2 (either ortho or para) excited by collisions with H2 (either ortho or
    para). The second argument should be one of the four H2 collision coefficient sets. */
double coolH2H2Polynomial(double T, const std::vector<double>& coefficients);

}; // namespace GloverAbel08

#endif // CORE_GLOVERABEL08_HPP
