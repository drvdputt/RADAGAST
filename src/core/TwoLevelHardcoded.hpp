#ifndef CORE_TWOLEVELHARDCODED_H_
#define CORE_TWOLEVELHARDCODED_H_

#include "NLevel.hpp"

/** A toy subclass of \c LevelDataProvider, for testing and simple demonstrations. */
class TwoLevelHardcoded : public LevelCoefficients
{
public:
	/** Constructs a \c LevelDataProvider object which will return the data to simulate a
	    toy model of CII 158 um. The data comes from
	    https://www.astro.umd.edu/~jph/N-level.pdf, bottom of page 4. */
	TwoLevelHardcoded();
	EMatrix cvv(const GasStruct& gas) const override;
};

#endif /* CORE_TWOLEVELHARDCODED_H_ */
