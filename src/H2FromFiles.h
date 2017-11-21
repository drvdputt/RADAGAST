#ifndef H2FROMFILES_H_
#define H2FROMFILES_H_

#include "LevelDataProvider.h"

#include <vector>

/** This class will implement the reading in and processing of all the files for the level model of
    H2. */
class H2FromFiles : public LevelDataProvider
{
public:
	H2FromFiles();
	void readData();
	int numLv() const override;
	EVector ev() const override;
	EVector gv() const override;
	EMatrix avv() const override;
	EMatrix extraAvv() const override;
	EMatrix cvv(double T, double ne, double np) const override;
	EVector sourcev(double T, double ne, double np) const override;
	EVector sinkv(double T, double n, double ne, double np) const override;

	enum class ElectronicState
	{
		X, B, Bprime, Cplus, Cminus, Dplus, Dminus
	};

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
		int e() const { return _e; }
	private:
		ElectronicState _eState;
		int _j, _v;
		double _e;
	};

private:
	std::vector<H2Level> _levelv;
};

#endif /* H2FROMFILES_H_ */
