#ifndef GASMODULE_GIT_SRC_H2FROMFILES_H_
#define GASMODULE_GIT_SRC_H2FROMFILES_H_

#include "CollisionData.h"
#include "LevelDataProvider.h"

#include <array>
#include <limits>
#include <map>
#include <vector>

/** This class will implement the reading in and processing of all the files for the level model
    of H2. */
class H2FromFiles : public LevelDataProvider
{
	//------------------------------//
	// CONSTRUCTION, READ-IN, SETUP //
	//------------------------------//
public:
	/** Creates a new \c H2FromFiles object, reads in all the data. Optional arguments:
	    upper limits for vibrational and rotational numbers. The data goes to 31 for J, and
	    to 14 for v. */
	H2FromFiles(int maxJ = 99, int maxV = 99);

	/** A clear way to index the electronic states of molecular hydrogen */
	enum class ElectronicState
	{
		// 1s 1Sigma+u
		X = 0,
		// 2p 1Sigma+u (Lyman)
		B,
		// 2p 1Pi+-u (C+ is Werner)
		Cplus,
		Cminus,
		// 3p 1Sigma+u
		Bprime,
		// 3p 1Pi+-u
		Dplus,
		Dminus
	};

private:

	/** Load the level energies. Comment in this file says 'by Evelyne Roueff'. There have
	    been some corrections, but the original data comes from Dabrowski (1984). */
	void readLevels();

	/** Load the data from Wolniewicz (1998) for X and abgrall (1994) for excited states.
	    See also MOLAT database. */
	void readTransProb();

	/** Load data from Lique (2015). */
	void readCollisions();

	/** Load dissociation cross sections from Gay et al. (2012). */
	void readDirectDissociation();

private:
	class H2Level
	{
	public:
		H2Level(ElectronicState eState, int j, int v, double e)
		                : _eState(eState), _j(j), _v(v), _e(e)
		{
		}
		ElectronicState eState() const { return _eState; }
		int j() const { return _j; }
		int v() const { return _v; }
		double e() const { return _e; }
		int g() const
		{
			// ortho = nuclear spin triplet / para = nuclear sping singlet
			bool ortho;
			bool oddJ = _j % 2;

			/* For these states, odd J -> ortho; even J -> para. */
			if (_eState == ElectronicState::X ||
			    _eState == ElectronicState::Cminus ||
			    _eState == ElectronicState::Dminus)
				ortho = oddJ;
			/* For the other states, it's the other way round. */
			else
				ortho = !oddJ;

			return ortho ? 3 * (2 * _j + 1) : 2 * _j + 1;
		}

	private:
		ElectronicState _eState;
		int _j, _v;
		double _e;
	};

	//-----------------------------------------------//
	// PUBLIC FUNCTIONS AKA THE OUTPUT OF THIS CLASS //
	//-----------------------------------------------//
public:
	size_t numLv() const override;

	/** Retrieve the index of a level with these quantum numbers. */
	size_t indexOutput(ElectronicState eState, int j, int v) const;

	/** The required overrides of the super class. */
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(double T, const EVector& speciesNv) const override;
	EVector sourcev(double T, const EVector& speciesNv) const override;
	EVector sinkv(double T, double n, const EVecto000r& speciesNv) const override;

	/** Functionality specific for H2. */

	/** Cross section for direct dissociation from @f$ X(J,v) @f$. */
	double directDissociationCrossSection(double nu, int j, int v);

private:
	/** Returns true if the given J and V are within the boundaries specified by the
	    user. */
	bool validJV(int J, int v) const;

	/** Adds the Cif and Cfi derived from the collision coefficient in qdata to the_cvv(i,
	    f) and the_cvv(f, i) respectively. For each transition in the CollisionData object,
	    q_if [cm3 s-1] is interpolated for the given temperature, and multiplied by the
	    given partner density to obtain the Cif [s-1]. Cfi is calculated using @f$ C_{fi} =
	    C_{if}\frac{g_i}{g_f}\exp(-h \nu_{if} / kT) @f$. */
	void addToCvv(EMatrix& the_cvv, const CollisionData& qdata, double T,
	              double nPartner) const;

	int _maxJ, _maxV;

	/* Contains the quantum numbers and energies of the levels. The output will be indexed
	   in the same way. */
	std::vector<H2Level> _levelv;
	size_t _numL{0};

	/* The dissociation energy of each electronic level. I assume the ionization threshold is 

	/* A map to help with converting (eState, j, v) quantum numbers to an index in the level
	   vector above. Function below shows how to add level. */
	std::map<std::array<int, 3>, int> _ejvToIndexm;
	inline void addLevel(ElectronicState eState, int j, int v, double e)
	{
		_levelv.emplace_back(eState, j, v, e);
		_ejvToIndexm.insert({{static_cast<int>(eState), j, v}, _levelv.size() - 1});
	}

	/** The transition coefficient matrix, indexed in the same way as _levelv, thanks to our
	    map approach. */
	EMatrix _avv;

	/** Collisions between H and H2. */
	// Cache species index of the collision partner
	int _inH;
	CollisionData _qH;

	// TODO: dissociating collisions
};

#endif /* GASMODULE_GIT_SRC_H2FROMFILES_H_ */
