#include "GloverAbel08.hpp"
#include "TemplatedUtils.hpp"

double GloverAbel08::coolOrthoH(double T)
{
	double T3 = T / 1000.;
	if (T < 100)
	{
		// equation 28
		return 5.09e-27 * std::sqrt(T3) * std::exp(-852.5 / T);
	}

	// equation 27
	double logT3 = std::log10(T3);
	double logCool = 1.;
	if (T < 1000)
	{
		logCool = TemplatedUtils::evaluatePolynomial(logT3, ortho_H_100to1000K_av);
	}
	else if (T < 6000)
	{
		logCool = TemplatedUtils::evaluatePolynomial(logT3, ortho_H_1000to6000K_av);
	}
	else
	{
		// For greater than 6000, use the result at 6000. This is also done in Cloudy
		// and in Gong et al. (2018).
		const double log6 = std::log10(6.);
		logCool = TemplatedUtils::evaluatePolynomial(log6, ortho_H_1000to6000K_av);
	}
	return std::pow(10., logCool);
}

double GloverAbel08::coolParaH(double T)
{
	double T3 = T / 1000.;
	if (T < 100)
	{
		// equation 29
		return 8.16e-26 * std::sqrt(T3) * std::exp(-509.85 / T);
	}

	// equation 27. Same as for ortho, but with different coefficients
	double logT3 = std::log10(T3);
	double logCool = 1.;
	if (T < 1000)
	{
		logCool = TemplatedUtils::evaluatePolynomial(logT3, para_H_100to1000K_av);
	}
	else if (T < 6000)
	{
		logCool = TemplatedUtils::evaluatePolynomial(logT3, para_H_1000to6000K_av);
	}
	else
	{
		// For greater than 6000, use the result at 6000. This is also done in Cloudy
		// and in Gong et al. (2018).
		const double log6 = std::log10(6.);
		logCool = TemplatedUtils::evaluatePolynomial(log6, para_H_1000to6000K_av);
	}
	return std::pow(10., logCool);
}

double GloverAbel08::coolParaPara(double T) { return coolH2H2Polynomial(T, para_para_av); }

double GloverAbel08::coolParaOrtho(double T) { return coolH2H2Polynomial(T, para_ortho_av); }

double GloverAbel08::coolOrthoPara(double T) { return coolH2H2Polynomial(T, ortho_para_av); }

double GloverAbel08::coolOrthoOrtho(double T) { return coolH2H2Polynomial(T, ortho_ortho_av); }

double GloverAbel08::coolH2H2Polynomial(double T, const std::vector<double>& coefficients)
{
	// Here I use a general function to evaluate one of the H2-H2 collision polynomials.
	// Otherwise, I would need to copy paste this temperature-dependent logic 4 times.
	if (T < 100)
		return 0.;

	double T3 = T / 1000.;
	double logT3 = std::log10(T3);
	double logCool = 1.;
	if (T < 6000)
	{
		// equation 31
		logCool = TemplatedUtils::evaluatePolynomial(logT3, coefficients);
	}
	else
	{
		const double log6 = std::log10(6.);
		logCool = TemplatedUtils::evaluatePolynomial(log6, coefficients);
	}
	return std::pow(10., logCool);
}
