#ifndef GASMODULE_GIT_SRC_H2LEVELS_H_
#define GASMODULE_GIT_SRC_H2LEVELS_H_

#include "NLevel.h"

class H2FromFiles;

class H2Levels: public NLevel
{
public:
	H2Levels(std::shared_ptr<const H2FromFiles> hff, const Array& frequencyv);
	~H2Levels();

	double dissociationRate(const Solution& s) const;
	double dissociationHeating(const Solution& s) const;

private:
	std::shared_ptr<const H2FromFiles> _hff;
};

#endif /* GASMODULE_GIT_SRC_H2LEVELS_H_ */
