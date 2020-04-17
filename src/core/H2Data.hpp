#ifndef CORE_H2FROMFILES_HPP
#define CORE_H2FROMFILES_HPP

#include "CollisionData.hpp"
#include "LevelCoefficients.hpp"
#include "Spectrum.hpp"
#include <array>
#include <limits>
#include <map>
#include <vector>

namespace GasModule
{
    /** This class implements the reading in and processing of all the files for the big H2
        model. */
    class H2Data : public LevelCoefficients
    {
        //------------------------------//
        // CONSTRUCTION, READ-IN, SETUP //
        //------------------------------//
    public:
        /** Create a new instance, reading in all the data. Optional arguments: upper limits for
            vibrational and rotational numbers. The data goes to 31 for J, and to 14 for v. */
        H2Data(int xMaxJ = 99, int xMaxV = 99, int eMaxJ = 99, int eMinV = 0, int eMaxV = 99);

        /** A clear way to index the electronic states of molecular hydrogen. The files from cloudy
            index these in the same order, from 0 to 6. To get this numerical index, it should be
            safe to simply do a static cast to an int. I also do that to make a map from the
            quantum numbers to the index of each level. */
        enum class ElectronicState {
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
        /** Load the level energies. Comment in this file says 'by Evelyne Roueff'. There have been
            some corrections, but the original data comes from Dabrowski (1984). */
        void readLevels();

        /** Loop over the loaded levels, and put all the energies in one vector. */
        EVector makeEv() const;

        /** Loop over the loaded levels, and put all the degeneracies in one vector. */
        EVector makeGv() const;

        /** Load the data from Wolniewicz (1998) (electric) and Pachucki and Komasa (2011)
            (magnetic) for X; Abgrall (1994) for excited states. See also MOLAT database. */
        void readTransProbs();

        /** Load data for the dissociation probabilities and the resulting kinetic energy for each
            level. */
        void readDissProbs();

        /** Load data from Lique (2015) for H, Lee (2008) for H2 (ortho and para), Gerlich (1990)
            for H+. */
        void readCollisions();

        /** Load dissociation cross sections from Gay et al. (2012). */
        void readDirectDissociation();

        /** Read energy levels and quantum numbers from the given file, registers them in the index
            map. */
        void readLevelFile(const std::string& repoFile, ElectronicState eState);

        /** Load the radiative transition rates between the levels from the given file. To be used
            after al levels have been read in. */
        void readTransProbFile(const std::string& repoFile, ElectronicState upperE, ElectronicState lowerE);

        /** Load the dissociation probabilities and kinetic energies from the given file. To be
            used after al levels have been read in. */
        void readDissProbFile(const std::string& repoFile, ElectronicState eState);

        // We don't have helium at this point, so leave it out. Don't forget to set numPartners
        // if you change this enum, as well as _gbarcoll down below.
        enum CollisionPartner {
            H0 = 0,
            // He,
            H2ORTHO,
            H2PARA,
            HPLUS
        };
        static const int _numPartners{4};

        /** Load collision data and store it in the given CollisionData object. To be used after al
            levels have been read in. */
        void readCollisionFile(const std::string& repoFile, CollisionPartner iPartner);

    public:
        /** This class contains the details about a single H2 level, and provides functions to
            calculate some simple derived properties of that level. */
        class H2Level
        {
        public:
            H2Level(ElectronicState eState, int j, int v, double e) : _eState(eState), _j(j), _v(v), _e(e)
            {
                bool oddJ = _j % 2;
                // For these states, odd J -> ortho; even J -> para. */
                if (_eState == ElectronicState::X || _eState == ElectronicState::Cminus
                    || _eState == ElectronicState::Dminus)
                    _ortho = oddJ;
                // For the other states, it's the other way round. This is explained in
                // the Cloudy H2 paper (Shaw et al. 2005).
                else
                    _ortho = !oddJ;
            }
            ElectronicState eState() const { return _eState; }
            int j() const { return _j; }
            int v() const { return _v; }
            double e() const { return _e; }
            int g() const
            {
                // ortho = nuclear spin triplet / para = nuclear spin singlet. Also,
                // remember that the states of a harmonic oscillator are not degenerate
                // in 1D, therefore v doesn't play a role for determining the
                // degeneracy.
                return _ortho ? 3 * (2 * _j + 1) : 2 * _j + 1;
            }
            bool ortho() const { return _ortho; }

        private:
            ElectronicState _eState;
            int _j, _v;
            double _e;
            bool _ortho;
        };

        /** Implement this inherited function to provide collision coefficients for the level
            transitions */
        EMatrix cvv(const CollisionParameters& cp) const override;

        /** Return the index of the level with the given quantum numbers. Returns -1 if level is
            not found. */
        int indexFind(ElectronicState eState, int j, int v) const;

        /** Return details of level at the given level index. */
        const H2Level& level(int index) const { return _levelv[index]; }

        /** Get the index pointing to the first (index-wise, not energy-wise) electronically
            excited level. All indices @f$ i < @f$ @c startOfExcitedIndices() correspond to levels
            of the ground state X, while @c startOfExcitedIndices() @f$ <= i <= @f$ @c numLv()
            point to levels of any of the electronically excited states. */
        size_t startOfExcitedIndices() const { return _startOfExcitedIndices; }

        /** The spontaneous dissociation probability per unit time from each level. [s-1] */
        const EVector& dissociationProbabilityv() const { return _dissProbv; }

        /** The average kinetic energy resulting from a spontaneous dissociation from each level.
            [erg] */
        const EVector& dissociationKineticEnergyv() const { return _dissKinEv; }

        /** Get the cross sections for direct dissociation from the given level. The Spectrum
            class is used for each cross section, so that the frequencies and cross sections are
            packaged together. */
        const std::vector<Spectrum>& directDissociationCrossSections(int index) const;

        /** Return a list of all the levels for which a direct dissociation cross section is
            available */
        const std::vector<int>& levelsWithCrossSectionv() const { return _levelsWithCrossSectionv; }

        /** Return the distribution of newly formed hydrogen into each level of the electronic
            ground state. Normalized to 1. For now, we use equation 19 from Draine and Bertoldi
            (1996), which does not depend on grain properties. Alternatively, the recipe from
            Takahashi (2001) could be used. */
        EVector formationDistribution() const;

    private:
        /** Returns true if the given J and V are within the boundaries specified by the user. */
        bool validJV(ElectronicState eState, int J, int v) const;

        /** Adds the Cif and Cfi derived from the collision coefficient in qdata to the_cvv(i, f)
            and the_cvv(f, i) respectively. For each transition in the CollisionData object, q_if
            [cm3 s-1] is interpolated for the given temperature, and multiplied by the given
            partner density to obtain the Cif [s-1]. Cfi is calculated using @f$ C_{fi} =
            C_{if}\frac{g_i}{g_f}\exp(-h \nu_{if} / kT) @f$. Any missing collisional data are
            approximated with the g-bar coefficients (which do not depend on temperature by the
            way). */
        void addToCvv(EMatrix& the_cvv, double T, CollisionPartner iPartner, double nPartner) const;

        /** Fill in a member containing the g-bar downward collision rates */
        void precalcGBarKvv();

        /** Uses the g-bar approximation (same as Cloudy, see Shaw et al. 2005, eq. 1 and Table 2)
            for the collision coefficients within X. This is a very rough approximation, but much
            better than just using 0 (personal communication with Peter van Hoof). The iPartner
            argument indicates which set of coefficients from Table 2 needs to be used, as a
            separate fit was made for each collision partner. See enum defined above. nPartner
            should be the density of this species. */
        void addGBarCvv(EMatrix& the_cvv, double kT, CollisionPartner iPartner, double nPartner) const;

        /** Returns the low-to-high collision coefficient Cfi, given the high to low coefficient
            Cif and the temperature kT. */
        double otherDirectionC(double Cif, int i, int f, double kT) const;

        // Settings
        bool _bB{false}, _bCplus{false}, _bCminus{false};
        int _groundMaxJ, _groundMaxV;
        int _excitedMaxJ, _excitedMinV, _excitedMaxV;

        // Contains the quantum numbers and energies of the levels
        std::vector<H2Level> _levelv;

        // The total number of levels (equivalent to _levelv.size())
        size_t _numL{0};
        size_t _startOfExcitedIndices{0};

        // A map to help with converting (eState, j, v) quantum numbers to an index in the level
        // vector above. Function below shows how adding a level works.
        std::map<std::array<int, 3>, int> _ejvToIndexm;
        inline void addLevel(ElectronicState eState, int j, int v, double e)
        {
            size_t newIndex = _levelv.size();
            _levelv.emplace_back(eState, j, v, e);
            _ejvToIndexm.insert({{static_cast<int>(eState), j, v}, newIndex});
        }

        // The transition coefficient matrix, indexed in the same way as _levelv
        EMatrix _avv;

        // Dissociation probability for each level. [s-1]
        EVector _dissProbv;
        // Kinetic energy following dissociation. [erg]
        EVector _dissKinEv;

        // Collision data for different partners, indexed using CollisionPartner enum.
        std::vector<CollisionData> _qdataPerPartner;

        // Keep track for which transition we have collisional data, for each partner individually.
        // If an entry _hadQdata[parter][i, f] is false, then this rate will be approximated using
        // gbar
        std::vector<EMatrix_bool> _hasQdata;

        // Coefficients (y0, a, b) to be used in the formula log(k) = y0 + a * max(sigma, 100)^b,
        // for each of the 5 possible collision partners defined in the enum above. (from Shaw et
        // al. 2005)
        const double _gbarcoll[_numPartners][3] = {{-9.9265, -0.1048, 0.456},
                                                   // {-8.281 , -0.1303 , 0.4931 },
                                                   {-10.0357, -0.0243, 0.67},
                                                   {-8.6213, -0.1004, 0.5291},
                                                   {-9.2719, -0.0001, 1.0391}};

        // Precalculate the downward g-bar rates. Indexed on [collision partner](upper, lower)
        std::vector<EMatrix> _gBarKvvPerPartner;

        // Cross sections for direct radiative dissociation from X. There can be multiple cross
        // sections for a level (first index), because there are different processes which can
        // occur per level (second index) (indicated by different 'nef' in the data file). */
        std::vector<std::vector<Spectrum>> _dissociationCrossSectionv;

        // Keep track of which levels have such cross sections, for easy looping
        std::vector<int> _levelsWithCrossSectionv;
    };
}
#endif  // CORE_H2FROMFILES_HPP
