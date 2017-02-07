#include <Eigen/Dense>
#include <vector>

class TwoLevel
{
public:
	TwoLevel();

	void doLevels(std::vector<double>& isrf, double T);

	std::vector<double> calculateEmission();

	std::vector<double> calculateOpacity();

	double n;		// Total density
	double n0, n1;	// Contribution to density of individual level populations

	// A matrix
	// formula for C matrix

};
