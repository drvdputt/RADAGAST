#ifndef _SRC_FREEFREE_H_
#define _SRC_FREEFREE_H_

#include "Table.h"

#include <string>
#include <vector>

class FreeFree
{
public:
	FreeFree();

private:
	void readData(std::string filename);

public:
	double gauntFactor(double temperature, double frequency) const;

private:
	// The ((natural!) log of the) data read from the file, indexed on (u, gamma^2)
	Table<2> _fileGauntFactorvv;
	double _loggamma2Min, _loguMin, _logStep, _loggamma2Max, _loguMax;

};

#endif /* _SRC_FREEFREE_H_ */
