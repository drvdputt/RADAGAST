#ifndef CORE_CHARGEDISTRIBUTION_HPP
#define CORE_CHARGEDISTRIBUTION_HPP

#include <functional>
#include <vector>

class ChargeDistribution
{
public:
    /** Default constructor for a Chargedistribution object. By default, a charge distribution with
        one entry is constructed, with charge 0 and density 1. This simply describes a population
        where all charges are 0. */
    ChargeDistribution() : _fz({1.}), _zmin{0} {}

    /** This function will calculate the detailed balance solution for the charge distribution,
        given two functions that produce the upward and downward charging rates for a certain z. A
        suggested charge range should also be given. The charge distribution will be cut off when
        the contributions become insignificant. The final zmin and zmax will lie within the given
        interval. */
    void calculateDetailedBalance(std::function<double(int z)> chargeUpRatef,
                                  std::function<double(int z)> chargeDownRatef, int zminGuess, int zmaxGuess);

    /** Minimum charge for which data is stored */
    int zmin() const { return _zmin; }

    /** Maximum charge for which data is stored */
    int zmax() const { return _zmin - 1 + _fz.size(); }

    /** Get the value of the charge distribution at the given z. Returns 0 if out of range. */
    double value(int z) const;

    /** Evaluate a function for every charge and sum the relative contributions. */
    double sumOverCharge(std::function<double(int z)> functionOfZ) const;

    /** Write the charge distribution to a two-column file (Z, fZ) */
    void plot(const std::string& file, const std::string& header) const;

private:
    std::vector<double> _fz;
    int _zmin;
};

#endif  // CORE_CHARGEDISTRIBUTION_HPP
