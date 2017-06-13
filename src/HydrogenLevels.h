#ifndef _SRC_HYDROGENLEVELS_H_
#define _SRC_HYDROGENLEVELS_H_

#include "NLevel.h"
#include "Table.h"

#include <vector>

class HydrogenLevels : public NLevel
{
public:
	HydrogenLevels();
	HydrogenLevels(const Array& frequencyv);

protected:
	int makeNLv() const override;
	Eigen::VectorXd makeEv() const override;
	Eigen::VectorXd makeGv() const override;
	Eigen::MatrixXd makeAvv() const override;
	Eigen::MatrixXd makeExtraAvv() const override;

	Eigen::MatrixXd prepareCollisionMatrix(double T, double electronDensity,
	                                       double protonDensity) const override;
	Array boundBoundContinuum(const Solution& s) const override;

private:
	/* Reads levels and their quantum numbers from CHIANTI, as well as the A coefficients
	 * between them. The listed levels are only those involved in spontaneous transitions, up to
	 * n = 5: 1s 2s 2p 2p 3s 3p 3p 3d 3d 4s 4p 4p 4d 4d 4f 4f 5s 5p 5p 5d 5d 5f 5f 5g 5g, and
	 * are also j-resolved (Hence the duplicates shown in this sentence). */
	void readData();

	/* Returns energy of a level read in from CHIANTI, given the principal (n) and angular
	 * momentum (l) numbers. Already averaged over different j. */
	double energyCHIANTI(int n, int l);

	/* Return the merged A coefficient using the values read in from the CHIANTI database. These
	 * are collapsed over the different j-values. */
	double einsteinACHIANTI(int ni, int li, int nf, int lf);
	/* With the initial level collapsed */
	double einsteinA(int ni, int li, int nf);
	/* With initial and final levels collapsed */
	double einsteinA(int ni, int nf);

	/* Return the the total electron collision strength (Upsilon). For the moment, there are
	 * only contributions from Anderson+2002 (J. Phys. B: At., Mol. Opt. Phys., 2002, 35, 1613).
	 * Need separate function for proton collision strength? */
	double eCollisionStrength(int ni, int li, int nf, int lf, double T);
	/* With the inital level collapsed */
	double eCollisionStrength(int ni, int nf, int lf, double T);
	/* With initial and final levels collapsed */
	double eCollisionStrength(int ni, int nf);

	typedef struct
	{
		int n, l, twoJplus1;
		double e;
	} levelInfo;

	const std::map<char, int> _lNumberm = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3}, {'G', 4}};

	/* Store the information about the levels in a vector. The index of a level in this vector
	is the number in the first column of the CHIANTI .elvlc file minus 1. */
	std::vector<levelInfo> _chiantiLevelv;

	/* To quickly find the level index for a set of quantum numbers, we use a map with fixed
	 * size arrays as keys {n, l, 2j+1} */
	std::map<std::array<int, 3>, int> _nljToChiantiIndexm;

	int _chiantiNumLvl{0};
	Table<2> _chiantiAvv;

	/* Indices like in Anderson 2000 table 1 */
	std::map<std::array<int, 2>, int> _nlToAndersonIndexm;

	// Temperature points, 8 of them in total, in electron volt
	Array _andersonTempv{{0.5, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0}};
	// Indexed on upper index, lower index, temperature index
	Table<3> _andersonCollStrvvv;

	bool _dataReady;

	/* -------some shorthands--- */

	/* Use this inline function for easy access */
	inline int indexCHIANTI(int n, int l, int twoJplus1)
	{
		return _chiantiLevelv[_nljToChiantiIndexm.at({n, l, twoJplus1})];
	}

	inline std::vector<int> twoJplus1range(int l)
	{
		return l > 0 ? std::vector<int>({2 * l, 2 * l + 2}) : std::vector<int>({2});
	}
};

#endif /* _SRC_HYDROGENLEVELS_H_ */
