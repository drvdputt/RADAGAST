#ifndef CORE_SPECIESINDEX_HPP
#define CORE_SPECIESINDEX_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"

#include <map>
#include <string>
#include <vector>

// class SpeciesIndex
// {
// private:
// 	/** This is a class with only static members. Constructor is private, so no instances can be
// 	    created. */
// 	SpeciesIndex();

// 	/** A map which maps species names to the indices used in the vector Currently available
// 	    names are "e-" for electrons "H+" for protons, "H" for atomic hydrogen and "H2" for
// 	    molecular hydrogen. */
// 	static const std::map<std::string, int> _indexMap;

// 	/** A function which generates the speciesIndex. */
// 	static std::map<std::string, int>
// 	createSpeciesIndexm(const std::vector<std::string>& namev);

// public:
// 	/** The index associated with a species name. Consult the source file to see what the
// 	    species names are. They will be used throughout the source code. */
// 	static int index(const std::string& s);

// 	/** Often used indices. */
// 	static int ine();
// 	static int inp();
// 	static int inH();
// 	static int inH2();

// 	/** The number of species. */
// 	static size_t size();

// 	/** Translates a list of species names (e.g. {'H', 'H+'}), plus corresponding
// 	    coefficients (e.g. {1, 1}) to a vector the size of the total number of species in
// 	    the simulation. The coefficients for the species not listed by the user are given
// 	    coefficient zero. In other words: translate the string for each species to unit
// 	    vector, and take a linear combination. */
// 	static EVector makeFullCoefficientv(const std::vector<std::string>& namev,
// 	                                    const Array& coefficientv);
// };

namespace SpeciesIndex
{
const std::vector<std::string> common4{"e-", "H+", "H", "H2"};
}

class NewSpeciesIndex
{
public:
	/** Create an empty species index. All index queries will return -1. */
	NewSpeciesIndex() = default;

	/** Create a species index, giving a list of names */
	NewSpeciesIndex(const std::vector<std::string>& namev);

	/** Add another species */
	void addSpecies(const std::string& name);

	/** Return -1 if not contained */
	int index(const std::string& name) const;

	/** Return vector, zero everywhere except at the index of the given species name, where
	    it is 1. */
	EVector unitVector(const std::string& name);

	/** Return number of species */
	size_t size() const { return _namev.size(); }

	/** Get species name for certain index */
	std::string name(int index) const { return _namev[index]; }

	/** Get all species names */
	std::vector<std::string> namev() const { return _namev; }

private:
	// Since it might be useful, let's just store both a map and a vector here. Vector gives
	// int --> string, map gives string --> int
	std::vector<std::string> _namev;
	std::map<std::string, int> _indexMap;
};

class SpeciesVector
{
public:
	/** Initialize with all elements zero */
	SpeciesVector(const NewSpeciesIndex& speciesIndex);

	/** Initialize and set elements */
	// SpeciesVector(const NewSpeciesIndex& speciesIndex, const EVector& nv);

	/** Things like the chemical network (which owns the SpeciesIndex) should be able to set
	    this directly. Note that this will crash if e-, H+, H or H2 are not prosent. */
	void setDensities(const EVector& nv) { _nv = nv; }
	void setNe(double ne) { _nv(_ine) = ne; }
	void setNp(double np) { _nv(_inp) = np; }
	void setNH(double nH) { _nv(_inH) = nH; }
	void setNH2(double nH2) { _nv(_inH2) = nH2; }

	double ni(int i) const { return i >= 0 ? _nv(i) : 0; }
	double ne() const { return ni(_ine); }
	double np() const { return ni(_inp); }
	double nH() const { return ni(_inH); }
	double nH2() const { return ni(_inH2); }
	double nSpecies(const std::string& name) const { return ni(_index->index(name)); }

	/** Useful for checking convergence etc. */
	EVector speciesNv() const { return _nv; }

	/** Useful for converting to other types */
	const double* data() const { return _nv.data(); }
	size_t size() const { return _nv.size(); }

private:
	const NewSpeciesIndex* _index; // storing as pointer is definitely more comfortable than
				       // as reference. Makes things copy-assignable.
	int _ine{-1}, _inp{-1}, _inH{-1}, _inH2{-1};
	EVector _nv;
};

#endif // CORE_SPECIESINDEX_HPP
