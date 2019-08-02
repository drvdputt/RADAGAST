#ifndef GASMODULE_GIT_SRC_CHARGEDISTRIBUTION_H_
#define GASMODULE_GIT_SRC_CHARGEDISTRIBUTION_H_

#include <functional>
#include <vector>

class ChargeDistribution
{
public:
	ChargeDistribution() : _fz({0.}), _zmin{0} {}

	/** This function will calculate the detailed balance solution for the charge
	    distribution, given two functions that produce the upward and downward charging
	    rates for a certain z. A suggested charge range should also be given. The charge
	    distribution will be cut off when the contributions become insignificant. The final
	    zmin and zmax will lie within the given interval. */
	void calculateDetailedBalance(std::function<double(int z)> chargeUpRatef,
	                              std::function<double(int z)> chargeDownRatef,
	                              int zminGuess, int zmaxGuess);

	int zmin() const { return _zmin; }
	int zmax() const { return _zmin - 1 + _fz.size(); }
	double value(int z) const;
	double sumOverCharge(std::function<double(int z)> functionOfZ) const;

private:
	std::vector<double> _fz;
	int _zmin;
};

#endif /* GASMODULE_GIT_SRC_CHARGEDISTRIBUTION_H_ */
