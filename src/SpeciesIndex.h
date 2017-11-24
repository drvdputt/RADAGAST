#ifndef GASMODULE_GIT_SRC_SPECIESINDEX_H_
#define GASMODULE_GIT_SRC_SPECIESINDEX_H_

#include <map>
#include <string>
#include <vector>

class SpeciesIndex
{
private:
	/** This is a class with only static members. Constructor is private, so no instances can be
	    created. */
	SpeciesIndex();
	
	/** A map which maps species names to the indices used in the vector Currently available
	    names are "e-" for electrons "H+" for protons, "H" for atomic hydrogen and "H2" for
	    molecular hydrogen. */
	static const std::map<std::string, int> indexMap;

	/** A function which generates the speciesIndex. */
	static std::map<std::string, int>
	createSpeciesIndexm(const std::vector<std::string>& namev);

public:
	/** The index associated with a species name. Consult the source file to see what the
	    species names are. They will be used throughout the source code. */
	static int index(const std::string& s);

	/** Often used indices. */
	static int ine();
	static int inp();
	

	/** The number of species. */
	static size_t size();
};

#endif /* GASMODULE_GIT_SRC_SPECIESINDEX_H_ */
