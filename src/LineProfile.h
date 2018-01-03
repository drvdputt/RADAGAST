#ifndef GASMODULE_GIT_SRC_LINEPROFILE_H
#define GASMODULE_GIT_SRC_LINEPROFILE_H

/** A set of tools to handle line profile in an efficient way. */
class LineProfile
{
public:
	LineProfile(double center, double sigma_gauss, double halfWidht_lorenz);
	// Evaluate the line profile at the given value (frequency) x
	double operator()(double nu) const;

private:
	double _center, _sigma_gauss;
	double _one_sqrt2sigma;
	double _a;
};
#endif /* GASMODULE_GIT_SRC_LINEPROFILE_H */
