#ifndef CORE_TWOLEVELHARDCODED_HPP
#define CORE_TWOLEVELHARDCODED_HPP

#include "LevelCoefficients.hpp"

/** A toy subclass, for testing and simple demonstrations. */
class TwoLevelHardcoded : public LevelCoefficients
{
public:
	/** Constructs an object which will contain the data to simulate a toy model of CII 158
	    um. The data comes from https://www.astro.umd.edu/~jph/N-level.pdf, bottom of page
	    4. */
	TwoLevelHardcoded();
	EMatrix cvv(const GasStruct& gas) const override;
};

#endif // CORE_TWOLEVELHARDCODED_HPP
