#ifndef CORE_TWOLEVELHARDCODED_H_
#define CORE_TWOLEVELHARDCODED_H_

#include "LevelDataProvider.hpp"

/** A toy subclass of \c LevelDataProvider, for testing and simple demonstrations. */
class TwoLevelHardcoded : public LevelDataProvider
{
public:
	/** Constructs a \c LevelDataProvider object which will return the data to simulate a
	    toy model of CII 158 um. The data comes from
	    https://www.astro.umd.edu/~jph/N-level.pdf, bottom of page 4. */
	TwoLevelHardcoded();

	/** Returns 2. */
	size_t numLv() const override;
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(const GasStruct& gas) const override;

private:
	EVector the_ev{EVector::Zero(2)};
	EVector the_gv{EVector::Zero(2)};
};

#endif /* CORE_TWOLEVELHARDCODED_H_ */
