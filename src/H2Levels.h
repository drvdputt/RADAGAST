#ifndef GASMODULE_GIT_SRC_H2LEVELS_H_
#define GASMODULE_GIT_SRC_H2LEVELS_H_

#include "NLevel.h"

class H2FromFiles;

class H2Levels: public NLevel
{
public:
	H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv);
	~H2Levels();

	/** The dissociation rate, both by direct photodissociation and the indirect Solomon process
	    derived from the level population solution. */
	double dissociationRate(const Solution& s, const Array& specificIntensityv) const;
	double dissociationHeating(const Solution& s) const;

private:
	std::shared_ptr<const H2FromFiles> _hff;
};

#endif /* GASMODULE_GIT_SRC_H2LEVELS_H_ */
