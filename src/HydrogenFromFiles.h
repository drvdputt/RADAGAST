#ifndef GASMODULE_GIT_SRC_HYDROGENFROMFILES_H_
#define GASMODULE_GIT_SRC_HYDROGENFROMFILES_H_

#include "LevelDataProvider.h"
#include "Table.h"

#include <map>
#include <vector>

class HydrogenFromFiles : public LevelDataProvider
{
	//------------------------------//
	// CONSTRUCTION, READ-IN, SETUP //
	//------------------------------//
public:
	HydrogenFromFiles(int resolvedUpTo = 5);

private:
	/* Reads levels and their quantum numbers from CHIANTI, as well as the A coefficients
	   between them. The listed levels are only those involved in spontaneous transitions, up to
	   n = 5: 1s 2s 2p 2p 3s 3p 3p 3d 3d 4s 4p 4p 4d 4d 4f 4f 5s 5p 5p 5d 5d 5f 5f 5g 5g, and
	   are also j-resolved (Hence the duplicates shown in this sentence). */
	void readData();

	/* Processes some of the data that was read to help with determining the final output. The
	   final number of levels, the parameters (n, l if resolved), and the order of the levels in
	   the output vectors and matrices are determined. All of these properties depend on the
	   number of resolved levels that was requested during construction. */
	void prepareForOutput();

	/* A struct-like class to store info about levels in a named way. When a quantum number has
	   the value -1, this means the level is collapsed over these numbers. To safeguard this
	   mechanism, I changed this to a class instead of a struct, to make sure the members stay
	   constant. */
	class HydrogenLevel
	{
	public:
		// Three constructors, for different degrees of collapsedness
		HydrogenLevel(int n, int l, int twoJplus1, double e)
		                : _n(n), _l(l), _twoJplus1(twoJplus1), _e(e)
		{
		}
		HydrogenLevel(int n, int l, double e) : _n(n), _l(l), _twoJplus1(-1), _e(e) {}
		HydrogenLevel(int n, double e) : _n(n), _l(-1), _twoJplus1(-1), _e(e) {}

		int n() const { return _n; }
		int l() const { return _l; }
		int twoJplus1() const { return _twoJplus1; }
		double e() const { return _e; }
		double g() const
		{
			if (_l < 0)
				return 2 * _n * _n;
			else if (_twoJplus1 < 0)
				return 4 * _l + 2;
			else
				return _twoJplus1;
		}

	private:
		int _n, _l, _twoJplus1; // quantum numbers
		double _e;              // energy
	};

	//-----------------------------------------------//
	// PUBLIC FUNCTIONS AKA THE OUTPUT OF THIS CLASS //
	//-----------------------------------------------//
	/* The output of these functions depends on the choice of the number of resolved levels.
	   They should be safe to call once prepareForOutput() has finished. The latter can only be
	   called after readData() has finished. These two calls should be made consecutively in
	   either the constructor or a setup routine. */
public:
	// FUNCTIONS RETURNING CONSTANT DATA //
	/* These functions convert the internally loaded data to usable matrices an vectors for
	   level calculations. These are usually called only once, during the setup of a typical
	   run. */

	/* Returns the number of levels */
	int numLv() const override;

	/* Returns a vector containing the energy of each level */
	Eigen::VectorXd ev() const override;

	/* Returns a vector containing the degeneracy of each level */
	Eigen::VectorXd gv() const override;

	/* Returns a matrix containing the Einstein A coefficients for all levels. Indexed on
	   (upper, lower), making it a lower triangle matrix. */
	Eigen::MatrixXd avv() const override;

	/* Returns a matrix containing any extra spontaneous decays between levels. This matrix can
	   be used to describe spontaneous decays that do NOT produce line radiation (for example
	   two-photon processes, which generate a continuum instead). */
	Eigen::MatrixXd extraAvv() const override;

	// FUNCTIONS RETURNING VARIABLE DATA //
	/* These functions provide coefficients that depend on external variables such as the
	   temperature. */

	/* Returns a matrix containing the collisional transition rates (already multiplied with the
	   partner density), for a given temperature and proton and electron densities. */
	Eigen::MatrixXd cvv(double T, double ne, double np) const override;

	/* Returns a vector containing the partial recombination rates into each level. */
	Eigen::VectorXd alphav(double T) const override;

	//-----------------------------------------//
	// FUNCTIONS THAT PROCESS THE READ-IN DATA //
	//-----------------------------------------//
	/* All of the functions below are private, and derive their output from 'scratch', meaning
	   directly from the data that was read-in during readData(). Once readData() has finished,
	   these should be safe to call. */
private:
	/* Returns energy of a level read in from CHIANTI, given the principal (n) and angular
	   momentum (l) numbers. Already averaged over different j. */
	double energy(int n, int l) const;
	/* Collapsed version */
	double energy(int n) const;

	/* Return the merged A coefficient using the values read in from the CHIANTI database. These
	   are collapsed over the different j-values. */
	double einsteinA(int ni, int li, int nf, int lf) const;
	/* With the initial level collapsed */
	double einsteinA(int ni, int nf, int lf) const;
	/* With initial and final levels collapsed */
	double einsteinA(int ni, int nf) const;

	/* Return the the total electron collision strength (Upsilon). For the moment, there are
	   only contributions from Anderson+2002 (J. Phys. B: At., Mol. Opt. Phys., 2002, 35, 1613).
	   Note that only queries for downward (in energy) transitions have the potential to return
	   a nonzero result.  Need separate function for proton collision strength? */
	double eCollisionStrength(int ni, int li, int nf, int lf, double T_eV) const;
	/* With the initial level collapsed */
	double eCollisionStrength(int ni, int nf, int lf, double T_eV) const;
	/* With initial and final levels collapsed */
	double eCollisionStrength(int ni, int nf, double T_eV) const;
	/* Works with LevelInfo objects, determining automatically what version to pick. Again, this
	   only works for downward transitions. */
	double eCollisionStrength(const HydrogenLevel& initial, const HydrogenLevel& final,
	                          double T_eV) const;

	/* A simple map for translating orbital letters into numbers */
	const std::map<char, int> _lNumberm = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3}, {'G', 4}};

	/* The total number of levels listed in the elvlc file from CHIANTI */
	int _chiantiNumLvl{0};

	/* Store the information about the levels in a vector. The index of a level in this vector
	  is the number in the first column of the CHIANTI .elvlc file minus 1. */
	std::vector<HydrogenLevel> _chiantiLevelv;

	/* To quickly find the level index as listed in the CHIANTI elvlc file for a set of quantum
	   numbers, we use a map with fixed size arrays as keys {n, l, 2j+1} */
	std::map<std::array<int, 3>, int> _nljToChiantiIndexm;

	/* With the following shorthand */
	inline int indexCHIANTI(int n, int l, int twoJplus1) const
	{
		return _nljToChiantiIndexm.at({n, l, twoJplus1});
	}

	/* The Einstein A coefficientes read in from the wgfa file from CHIANTI */
	Eigen::MatrixXd _chiantiAvv;

	/* Indices like in Anderson 2000 table 1 */
	std::map<std::array<int, 2>, int> _nlToAndersonIndexm;

	/* Temperature points, 8 of them in total, in electron volt */
	Array _andersonTempv{{0.5, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0}};

	/* Indexed on upper index, lower index, temperature index. Only downward collisions are
	   listed, so only store those, using a map which uses a pair of ints representing the upper
	   and lower levels. */
	std::map<std::array<int, 2>, Array> _andersonUpsilonvm;

	//----------------------//
	// BASED ON USER CHOICE //
	//----------------------//
	/* Below are some variables that help with providing the output. These variables and the
	   final output both depend on the configuration chosen by the user, such as the n up to
	   which the levels are resolved, declared below. */
private:
	int _resolvedUpTo{5};

	/* Number of levels outputted by this LevelDataProvider */
	int _numL{0};

	/* Contains the quantum numbers of the levels actually used for the output. l = -1 means
	   that the level is collapsed. The energy levels, A-coefficients, ... outputted will be
	   indexed the same way as in this vector. Changing the way the vector below is filled will
	   allow for some customization of the ordering of the levels. */
	std::vector<HydrogenLevel> _levelOrdering;
};

#endif /* GASMODULE_GIT_SRC_HYDROGENFROMFILES_H_ */
