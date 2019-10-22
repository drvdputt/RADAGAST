#ifndef CORE_HYDROGENHARDCODED_HPP
#define CORE_HYDROGENHARDCODED_HPP

#include "HData.hpp"

/** This provides data for a hydrogen level model in a hardcoded way. It is mainly for educational
    and compararive purposes. I used it to check if the implementation of \c HydrogenFromFiles
    correctly processed the j-resolved coefficients it read in, for example. */
class HydrogenHardcoded : public HData
{
public:
	HydrogenHardcoded();

	int nMax() const override { return 5; }

	size_t index(int n, int l) const override;

	std::array<size_t, 2> twoPhotonIndices() const override;

	EMatrix cvv(const GasStruct& gas) const override;

private:
	EVector the_ev;
	EVector the_gv;
};

#endif // CORE_HYDROGENHARDCODED_HPP
