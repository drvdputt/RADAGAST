#include "ChargeDistribution.hpp"
#include "Error.hpp"

double ChargeDistribution::value(int z) const
{
	if (z < _zmin || z > zmax())
		return 0;
	else
		return _fz[z - _zmin];
}

double ChargeDistribution::sumOverCharge(std::function<double(int z)> functionOfZ) const
{
	double sum = 0;
	for (int i = 0; i < _fz.size(); i++)
	{
		if (!std::isfinite(_fz[i]))
			Error::runtime("nan in charge distribution");
		sum += _fz[i] * functionOfZ(_zmin + i);
	}
	return sum;
}

void ChargeDistribution::calculateDetailedBalance(std::function<double(int z)> chargeUpRatef,
                                                  std::function<double(int z)> chargeDownRatef,
                                                  int zminGuess, int zmaxGuess)
{
	_fz.resize(zmaxGuess - zminGuess + 1, -1);

	/* We will cut off the distribution at some point past the maximum (in either the
	   positive or the negative direction). For a Gaussian distribution, at half the height
	   we have already ~80% of the population. So a 10th of the height should definitely be
	   sufficient.*/
	double cutOffFactor = 1.e-1;
	int trimLow = 0;
	int trimHigh = _fz.size() - 1;
	int centerZ = 0;

	// Find the central peak using a binary search
	int lowerbound = zminGuess;
	int upperbound = zmaxGuess;
	int currentZ = floor((zminGuess + zmaxGuess) / 2.);
	while (currentZ != lowerbound)
	{
		// Compare the up and down rates for Z and Z+1 respectively
		// goal: up * f(z) == down * f(z+1)
		// factor up/down == f(z+1)/f(z)
		// --> factor > 1 means upward slope and vice versa
		double factor = chargeUpRatef(currentZ) / chargeDownRatef(currentZ + 1);
		if (factor > 1)
		{
			// Upward slope --> the maximum is more to the right
			lowerbound = currentZ;
		}
		else
		{
			// Downward slope --> the maximum is more to the left
			upperbound = currentZ;
		}
		// Move the cursor to the center of the new bounds
		currentZ = floor((lowerbound + upperbound) / 2.);
	}
	// The result of the binary search
	centerZ = currentZ;

	// Apply detailed balance equation: start at one ... The equation which must be
	// satisfied is f(Z) * upRate(Z) = f(Z+1) * downRate(Z+1)
	_fz[centerZ - zminGuess] = 1;

	// ... for Z > centerZ
	for (int z = centerZ + 1; z <= zmaxGuess; z++)
	{
		int index = z - zminGuess;
		// f(z) =  f(z-1) * up(z-1) / down(z)
		double up = chargeUpRatef(z - 1);
		double down = chargeDownRatef(z);
		_fz[index] = down > 0 ? _fz[index - 1] * up / down : 0;

		if (!std::isfinite(_fz[index]))
			Error::runtime("invalid value in charge distribution");

		if (_fz[index] < cutOffFactor)
		{
			trimHigh = index;
			break;
		}
	}

	// ... for Z < centerZ
	for (int z = centerZ - 1; z >= zminGuess; z--)
	{
		int index = z - zminGuess;
		// f(z) = f(z+1) * down(z+1) / up(z)
		double up = chargeUpRatef(z);
		double down = chargeDownRatef(z + 1);
		_fz[index] = up > 0 ? _fz[index + 1] * down / up : 0;

		if (!std::isfinite(_fz[index]))
			Error::runtime("invalid value in charge distribution");

		if (_fz[index] < cutOffFactor)
		{
			trimLow = index;
			break;
		}
	}

	/* Apply the cutoffs (this might be done more elegantly, i.e. while avoiding the
	   allocation of all the memory for the initial fZ buffer). */

	// Trim the top first! Makes things (i.e. dealing with the indices) much easier.

	_fz.erase(_fz.begin() + trimHigh + 1, _fz.end());
	_fz.erase(_fz.begin(), _fz.begin() + trimLow);
	_zmin = zminGuess + trimLow;

	// Normalize
	double sum = 0.;
	for (double d : _fz)
		sum += d;
	for (double& d : _fz)
		d /= sum;
}
