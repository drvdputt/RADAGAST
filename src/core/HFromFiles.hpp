#ifndef CORE_HYDROGENFROMFILES_HPP
#define CORE_HYDROGENFROMFILES_HPP

#include "Array.hpp"
#include "CollisionData.hpp"
#include "HData.hpp"
#include <map>
#include <vector>

namespace GasModule
{
    /** This class implements the reading in and processing of all the files for the 5-level model
        of H. */
    class HFromFiles : public HData
    {
        //------------------------------//
        // CONSTRUCTION, READ-IN, SETUP //
        //------------------------------//
    public:
        /** Creates a new instance, reading in all the data. Optional argument: the number levels
            up to which l is taken into account. */
        HFromFiles(int resolvedUpTo = 5);

    private:
        /** Reads levels and their quantum numbers from CHIANTI, as well as the A coefficients
            between them. The listed levels are only those involved in spontaneous transitions, up
            to n = 5: 1s 2s 2p 2p 3s 3p 3p 3d 3d 4s 4p 4p 4d 4d 4f 4f 5s 5p 5p 5d 5d 5f 5f 5g 5g,
            and are also j-resolved (Hence the duplicates shown in this sentence). */
        void readData();

        /** Processes some of the data that was read to help with determining the final
            coefficients. The final number of levels, the parameters (n, l if resolved), and the
            order of the levels are determined. All of these properties depend on the number of
            resolved levels that was requested during construction. */
        void prepareForOutput();

        EVector makeEv() const;
        EVector makeGv() const;
        EMatrix makeAvv() const;
        EMatrix makeExtraAvv() const;

        /** A struct-like class to store info about levels in a named way. When a quantum number
            has the value -1, this means the level is collapsed over these numbers. To safeguard
            this mechanism, I changed this to a class instead of a struct, to make sure the members
            stay constant. */
        class HydrogenLevel
        {
        public:
            /** Constructor for an nlJ-resolved level. Takes \f$n\f$ and \f$l\f$ as integers, and
                \f$2j + 1\f$ instead of \f$j\f$ so an int can be used. */
            HydrogenLevel(int n, int l, int twoJplus1, double e) : _n(n), _l(l), _twoJplus1(twoJplus1), _e(e) {}
            /** Constructor for an nl-resolved level. Puts -1 at the place of \f$2j+1\f$. */
            HydrogenLevel(int n, int l, double e) : _n(n), _l(l), _twoJplus1(-1), _e(e) {}

            HydrogenLevel(int n, double e) : _n(n), _l(-1), _twoJplus1(-1), _e(e) {}

            bool nljResolved() const { return _twoJplus1 >= 0; }
            bool nlResolved() const { return _l >= 0 && _twoJplus1 < 0; }
            bool nCollapsed() const { return _l < 0; }
            int n() const { return _n; }
            int l() const { return _l; }
            int twoJplus1() const { return _twoJplus1; }
            double e() const { return _e; }

            /** Degeneracy of the level. Double because ratios of different g are often used. */
            double g() const
            {
                if (nCollapsed())
                    return 2 * _n * _n;
                else if (nlResolved())
                    return 4 * _l + 2;
                else
                    return _twoJplus1;
            }

        private:
            int _n, _l, _twoJplus1;  // quantum numbers
            double _e;               // energy
        };

    public:
        /** Calculate the collision rates. Contributions: 1. n-changing rates based on the
            electron collision strength data Anderson+2002 (J. Phys. B: At., Mol. Opt. Phys.,
            2002, 35, 1613). 2. l-changing collision rate coefficients as described by Pengelley
            \& Seaton (1964) */
        EMatrix cvv(const CollisionParameters& cp) const override;

        int nMax() const override { return 5; }

        /** Gives the index of the (n, l) energy level which is used in the output functions. To
            get the energy of the n=5, l=2 level for example, one can call i = indexOutput(5, 2).
            The energy you need will be the i'th element of the output vector produced by the
            ev()-function (which the client should have cached, since this function is "slow"). If
            the n'th level is collapsed, the given l will be ignored. */
        size_t index(int n, int l = 0) const override;

        /** Return a pair of indices indication the upper and lower level of the two-photon
            transition (2s ans 1s respectively). When the upper level is collapsed, the index
            corresponding to n=2 will be given. In this case, the extra transition rate at these
            indices, extraAvv(n=2, n=1), equals 1/4 of the two-photon transition rate, which is the
            same as assuming that 1/4 of the n=2 atoms is in a 2s state. Therefore, in case the n=2
            level is collapsed, the two photon continuum emission should be rescaled with the same
            factor. This will probably give very inaccurate results, but its the best you can do
            with a collapsed n=2 level. */
        std::array<size_t, 2> twoPhotonIndices() const override;

    private:
        /** Calculates the l-changing collision rate coefficients as described by Pengelley \&
            Seaton (1964) Since l can only change by 1, the calculation starts at either end (li=0
            or li=n-1) of the l-range, and calculates the q_nli->nlf based on the q_n(li-1),n(lf-1)
            or q_n(li+1),n(lf-1). The results are returned as a tridiagonal matrix indexed on
            (li,lf), of dimension n x n. Used only by cvv, and after prepareForOutput, and hence
            declared private here. [cm3 s-1] */
        EMatrix PS64CollisionRateCoeff(int n, double T, double ne) const;

        //---------------------------------------------//
        // FUNCTIONS DEALING WITH COLLAPSING OF LEVELS //
        //---------------------------------------------//

        /** Returns energy of a level read in from CHIANTI, given the principal (n) and angular
            momentum (l) numbers. Already averaged over different j. [erg] */
        double energy(int n, int l) const;
        /** Collapsed version */
        double energy(int n) const;

        // Chianti data
        // ------------

        // The levels, stored as they were read in (j-resolved). The index of a level in this
        // vector is the number in the first column of the CHIANTI .elvlc file minus 1.
        std::vector<HydrogenLevel> _chiantiLevelv;

        // The Einstein A coefficients read in from the wgfa file from CHIANTI
        EMatrix _chiantiAvv;

        // Map from quantum numbers to level index as listed in the CHIANTI elvlc file. Uses fixed
        // size arrays as keys {n, l, 2j+1}.
        std::map<std::array<int, 3>, int> _nljToChiantiIndexm;
        inline int indexCHIANTI(int n, int l, int twoJplus1) const { return _nljToChiantiIndexm.at({n, l, twoJplus1}); }

        // Anderson collision data
        // -----------------------

        // Correspondence between the in the Anderson data file (minus 1), and n,l. Minus 1
        // because the indices in the file start from 1.
        std::vector<std::array<int,2>> _andersonIndexm1ToNLv;

        // Store the data using this class
        CollisionData _qdataAnderson;

        // Misc
        // ----

        // Translation from orbital letters into numbers
        const std::map<char, int> _lNumberm = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3}, {'G', 4}};

        // Total spontaneous decay rate of each level. Needed for the l-changing collision formula
        // of PS64.
        EVector _totalAv;

        // Number of levels
        size_t _numL{0};

        // Highest n for which l is resolved
        int _resolvedUpTo{5};

        // Quantum numbers of the levels. l = -1 means that the level is collapsed. The energy
        // levels, A-coefficients are stored in the same order as this vector.
        std::vector<HydrogenLevel> _levelOrdering;

        // Map from {n, l} to the index in _levelOrdering.
        std::map<std::array<int, 2>, size_t> _nlToOutputIndexm;
    };
}
#endif  // CORE_HYDROGENFROMFILES_HPP
