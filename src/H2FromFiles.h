#ifndef GASMODULE_GIT_SRC_H2FROMFILES_H_
#define GASMODULE_GIT_SRC_H2FROMFILES_H_

#include "CollisionData.h"
#include "LevelDataProvider.h"
#include "Spectrum.h"

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

	/** A clear way to index the electronic states of molecular hydrogen. The files from
	    cloudy index these in the same order, from 0 to 6. To get this numerical index, it
	    should be safe to simply do a static cast to an int. I also do that to make a map
	    from the quantum numbers to the index of each level. */
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

	/** Load the data from Wolniewicz (1998) (electric) and Pachucki and Komasa (2011)
	    (magnetic) for X; Abgrall (1994) for excited states. See also MOLAT database. */
	void readTransProbs();

	/** Load data for the dissociation probabilities and the resulting kinetic energy for
	    each level. */
	void readDissProbs();

	/** Load data from Lique (2015), Lee (2008). */
	void readCollisions();

	/** Load dissociation cross sections from Gay et al. (2012). */
	void readDirectDissociation();

	void readLevelFile(const std::string& repoFile, ElectronicState eState);

	void readTransProbFile(const std::string& repoFile, ElectronicState upperE,
	                       ElectronicState lowerE);

	void readDissProbFile(const std::string& repoFile, ElectronicState eState);

	CollisionData readCollisionFile(const std::string& repoFile) const;

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
			bool ortho;
			bool oddJ = _j % 2;

			/* For these states, odd J -> ortho; even J -> para. */
			if (_eState == ElectronicState::X ||
			    _eState == ElectronicState::Cminus ||
			    _eState == ElectronicState::Dminus)
				ortho = oddJ;
			/* For the other states, it's the other way round. This is explained in
			   the Cloudy H2 paper (Shaw et al. 2005). */
			else
				ortho = !oddJ;

			// ortho = nuclear spin triplet / para = nuclear spin singlet. Also,
			// remember that the states of a harmonic oscillator are not degenerate
			// in 1D, therefore v doesn't play a role for determining the
			// degeneracy.
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
	/** @name The required overrides of the super class. */
	/**@{*/
	size_t numLv() const override;
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(const GasStruct& gas) const override;
	/**@}*/

	/** @name Functionality specific for H2. */
	/**@{*/

	/** Retrieve the index of a level with these quantum numbers. */
	size_t indexOutput(ElectronicState eState, int j, int v) const;

	/** Get a list of all the indices that correspond levels of the electronic ground (X)
	    state. */
	std::vector<size_t> indicesX() const { return _indicesX; }

	/** Get a list of all the indices that correspond levels of any of the electronic
	    excited states. */
	std::vector<size_t> indicesExcited() const { return _indicesExcited; }

	EVector dissociationProbabilityv() const { return _dissProbv; }
	EVector dissociationKineticEnergyv() const { return _dissKinEv; }

	/** Cross section for direct dissociation from @f$ X(J,v) @f$. */
	double directDissociationCrossSection(double nu, int j, int v) const;

	/** Cross section for direct dissociation from level with index @c index. */
	double directDissociationCrossSection(double nu, size_t index) const;

	/** Get all the cross sections as a vector of Spectrum objects. This is handy if you
	    want to integrate over these cross sections efficiently, since the Spectrum object
	    contains the minimum and maximum frequency. */
	std::vector<Spectrum> directDissociationCrossSections(size_t index) const;
	/**@}*/
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

	/* The total number of levels (equivalent to _levelv.size()) */
	size_t _numL{0};

	/* List of levels that belong to the electronic ground state, and any of the excited
	   states, respectively. */
	std::vector<size_t> _indicesX;
	std::vector<size_t> _indicesExcited;

	/* A map to help with converting (eState, j, v) quantum numbers to an index in the level
	   vector above. Function below shows how adding a level works. */
	std::map<std::array<int, 3>, int> _ejvToIndexm;
	inline void addLevel(ElectronicState eState, int j, int v, double e)
	{
		size_t newIndex = _levelv.size();
		_levelv.emplace_back(eState, j, v, e);
		_ejvToIndexm.insert({{static_cast<int>(eState), j, v}, newIndex});
		if (eState == ElectronicState::X)
			_indicesX.emplace_back(newIndex);
		else
			_indicesExcited.emplace_back(newIndex);
	}

	/** The transition coefficient matrix, indexed in the same way as _levelv, thanks to our
	    map approach. */
	EMatrix _avv;

	/** Collisions between H and H2. */
	CollisionData _qH;

	/** Collisions between H2 and H2. */
	CollisionData _qH2ortho, _qH2para;

	/** Dissociation probability for each level. [s-1] */
	EVector _dissProbv;
	/** Kinetic energy following dissociation. [erg] */
	EVector _dissKinEv;

	// Caches species index of the collision partners
	int _inH;
	int _inH2;

	/** Cross sections for direct radiative dissociation from X. Indexed on (level index,
	    process). There can be multiple cross sections for a level, because there are
	    different processes which can occur per level (indicated by different 'nef' (final
	    electronic level) in the data file. These will be indexed arbitrarily, as the nature
	    of the process shouldn't matter. */
	std::vector<std::vector<Spectrum>> _dissociationCrossSectionv;
};

#endif /* GASMODULE_GIT_SRC_H2FROMFILES_H_ */
