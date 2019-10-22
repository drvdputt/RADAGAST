#ifndef CORE_SPECIESINDEX_HPP
#define CORE_SPECIESINDEX_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

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
	static const std::map<std::string, int> _indexMap;

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
	static int inH();
	static int inH2();

	/** The number of species. */
	static size_t size();

	/** Translates a list of species names (e.g. {'H', 'H+'}), plus corresponding
	    coefficients (e.g. {1, 1}) to a vector the size of the total number of species in
	    the simulation. The coefficients for the species not listed by the user are given
	    coefficient zero. In other words: translate the string for each species to unit
	    vector, and take a linear combination. */
	static EVector makeFullCoefficientv(const std::vector<std::string>& namev,
	                                    const Array& coefficientv);
};

#endif // CORE_SPECIESINDEX_HPP
