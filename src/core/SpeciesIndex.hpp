#ifndef CORE_SPECIESINDEX_HPP
#define CORE_SPECIESINDEX_HPP

#include "Array.hpp"
#include "EigenAliases.hpp"
#include <map>
#include <string>
#include <vector>

class SpeciesIndex
{
public:
    /** Create an empty species index. All index queries will return -1. */
    SpeciesIndex() = default;

    /** Create a species index, giving a list of names */
    SpeciesIndex(const std::vector<std::string>& namev);

    /** Add another species */
    void addSpecies(const std::string& name);

    /** Return -1 if not contained */
    int index(const std::string& name) const;

    /** Sum of unit vectors for the given names, with the given weights. Useful to set
	    multiple densities, or to represent a reaction. */
    EVector linearCombination(const std::vector<std::string>& namev, const Array& weightv) const;

    /** Return number of species */
    size_t size() const { return _namev.size(); }

    /** Get species name for certain index */
    const std::string& name(int index) const { return _namev[index]; }

    /** Get all species names */
    const std::vector<std::string>& namev() const { return _namev; }

private:
    // Since it might be useful, let's just store both a map and a vector here. Vector gives
    // int --> string, map gives string --> int
    std::vector<std::string> _namev;
    std::map<std::string, int> _indexMap;

public:
    static const std::vector<std::string> e_p_H_H2;
};

class SpeciesVector
{
public:
    /** Initialize with all elements zero */
    SpeciesVector(const SpeciesIndex* speciesIndex);

    /** Things like the chemical network (which owns a SpeciesIndex) should be able to set
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

    const EVector& speciesNv() const { return _nv; }

    const double* data() const { return _nv.data(); }
    size_t size() const { return _nv.size(); }

    /** Overload the output stream operator, for easy printing and debugging. */
    friend std::ostream& operator<<(std::ostream& os, const SpeciesVector& sv);

private:
    const SpeciesIndex* _index;
    int _ine{-1}, _inp{-1}, _inH{-1}, _inH2{-1};
    EVector _nv;
};

#endif  // CORE_SPECIESINDEX_HPP
