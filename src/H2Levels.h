#ifndef GASMODULE_GIT_SRC_H2LEVELS_H_
#define GASMODULE_GIT_SRC_H2LEVELS_H_

#include "NLevel.h"

class H2FromFiles;

class H2Levels : public NLevel
{
public:
	H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv);
	~H2Levels();

	/** An override which is better suited for H2. It scales as a*n^2, where a is the number
	    of iterations, instead of n^3, apparently. I might be interesting to see this with
	    my own eyes. */
	EVector solveRateEquations(double n, const EMatrix& BPvv, const EMatrix& Cvv,
	                           const EVector& sourcev, const EVector& sinkv,
	                           int chooseConsvEq) const override;

	/** The dissociation rate, both by direct photodissociation and the indirect Solomon
	    process derived from the level population solution. */
	double dissociationRate(const Solution& s, const Array& specificIntensityv) const;
	double dissociationHeating(const Solution& s) const;

private:
	std::shared_ptr<const H2FromFiles> _hff;
};

#endif /* GASMODULE_GIT_SRC_H2LEVELS_H_ */
