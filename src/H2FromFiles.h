#ifndef H2FROMFILES_H_
#define H2FROMFILES_H_

#include "LevelDataProvider.h"

#include <map>
#include <vector>

/** This class will implement the reading in and processing of all the files for the level model of
    H2. */
class H2FromFiles : public LevelDataProvider
{
	//------------------------------//
	// CONSTRUCTION, READ-IN, SETUP //
	//------------------------------//
public:
	/** Creates a new \c H2FromFiles object, reads in all the data. */
	H2FromFiles();

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
	void readData();

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
			if (_eState == ElectronicState::X || _eState == ElectronicState::Cminus ||
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
	EVector ev() const override;

	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(double T, double ne, double np) const override;
	EVector sourcev(double T, double ne, double np) const override;
	EVector sinkv(double T, double n, double ne, double np) const override;

private:
	/* Contains the quantum numbers and energies of the levels. The output will be indexed in
	   the same way. */
	std::vector<H2Level> _levelv;
	size_t _numL{0};

	/* A map to help with converting (eState, j, v) quantum numbers to an index in the level
	   vector above. Function below shows how to add level. */
	std::map<std::array<int, 3>, int> _evjToIndex;
	inline void addLevel(ElectronicState eState, int j, int v, double e)
	{
		_levelv.emplace_back(eState, j, v, e);
		_evjToIndex.insert({{static_cast<int>(eState), j, v}, _levelv.size() - 1});
	}

	/* The transition coefficient matrix, indexed in the same way as _levelv, thanks to our map
	   approach. */
	EMatrix _avv;
};

#endif /* H2FROMFILES_H_ */
