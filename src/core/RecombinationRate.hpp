#ifndef GASMODULE_GIT_SRC_RECOMBINATIONRATE_H_
#define GASMODULE_GIT_SRC_RECOMBINATIONRATE_H_

#include "Array.hpp"
#include "Table.hpp"

#include <array>
#include <map>
#include <vector>

class RecombinationRate
{
public:
	RecombinationRate() = default;
	virtual ~RecombinationRate() = default;
	virtual double alpha(int n, int l, double T) const = 0;
};

class HydrogenADF48 : public RecombinationRate
{
public:
	HydrogenADF48();
	~HydrogenADF48();

private:
	void readADF48File(const std::string& path);

public:
	double alpha(int n, int l, double T) const override;

private:
	std::map<std::array<int, 2>, size_t> _nlToIndexm;
	Array _temperaturev;
	std::vector<std::vector<double>> _alphavv;
};
#endif /* GASMODULE_GIT_SRC_RECOMBINATIONRATE_H_ */
