#include "HFromFiles.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "Ionization.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"

using namespace std;

namespace
{
    const int cNMAX = 5;
    const string chiantiBaseName = "/dat/CHIANTI_8.0.6_data/h/h_1/h_1";

    inline vector<int> twoJplus1range(int l)
    {
        // If l > 0, then j can be either l + 0.5 or l - 0.5. If l==0, j is always 1/2 and 2j+1
        // can only be 2
        return l > 0 ? vector<int>({2 * l, 2 * l + 2}) : vector<int>({2});
    }
}

namespace GasModule
{
    HFromFiles::HydrogenLevel::HydrogenLevel(int n, int l, int twoJplus1, double e)
        : _n(n), _l(l), _twoJplus1(twoJplus1), _e(e)
    {}
    HFromFiles::HydrogenLevel::HydrogenLevel(int n, int l, double e) : _n(n), _l(l), _twoJplus1(-1), _e(e) {}
    HFromFiles::HydrogenLevel::HydrogenLevel(int n, double e) : _n(n), _l(-1), _twoJplus1(-1), _e(e) {}

    double HFromFiles::HydrogenLevel::g() const
    {
        if (_l == -1 && _twoJplus1 == -1)
        {
            return 2 * _n * _n;
        }
        else if (_l != -1)
        {
            if (_twoJplus1 == -1)
                return 4 * _l + 2;
            else
                return _twoJplus1;
        }
        else
        {
            Error::runtime("Wrong combination of n, l, and 2j+1 for H level");
            return -1;
        }
    }

    HFromFiles::HFromFiles(int resolvedUpTo) : _resolvedUpTo(resolvedUpTo)
    {
        if (_resolvedUpTo > cNMAX) Error::rangeCheck<int>("Number of resolved levels", _resolvedUpTo, 2, cNMAX);

        // read and process data in the correct order (processing of transition probabilities
        // and collision strengths might depend on the way the levels were set up)
        readLevels();
        readTransProbs();
        readCollisions();

        // Set required members of parent class
        setConstants(makeEv(), makeGv(), makeAvv(), makeExtraAvv());
    }

    void HFromFiles::readLevels()
    {
        // Map from quantum numbers to level index as listed in the CHIANTI elvlc file. Uses
        // fixed size arrays as keys {n, l, 2j+1}. We will use it at the end of this function to
        // process the nlj levels into nl levels.
        std::map<std::array<int, 3>, int> nljToChiantiIndexm;

        // Translation from orbital letters into numbers
        const std::map<char, int> lNumberm = {{'S', 0}, {'P', 1}, {'D', 2}, {'F', 3}, {'G', 4}};

        ifstream elvlc = IOTools::ifstreamRepoFile(chiantiBaseName + ".elvlc");
        string line;

        // Start with the first line
        getline(elvlc, line);
        while (line.compare(1, 2, "-1"))
        {
            // Read the different parts of the line
            int lvIndex, twoSplus1;
            string config;
            char lSymbol;
            double j, observedEnergy, theoreticalEnergy;
            istringstream(line) >> lvIndex >> config >> twoSplus1 >> lSymbol >> j >> observedEnergy
                >> theoreticalEnergy;

            // Get the first number from the config string
            int n;
            istringstream(config) >> n;

            // Translate the angular momentum letter
            int l = lNumberm.at(lSymbol);

            // Store 2j+1 (static cast to make it clear that j is not integer)
            int twoJplus1 = static_cast<int>(2 * j + 1);

            // Convert the energy from cm-1 to erg
            double e = observedEnergy * Constant::LIGHT * Constant::PLANCK;

            // The level indices in the data structures will go from 0 to number of levels minus
            // one. The level indices in the file go from 1 to the number of levels, while those
            // used for the vector and the map will go from 0 to numL - 1.
            _chiantiLevelv.emplace_back(n, l, twoJplus1, e);
            nljToChiantiIndexm.insert({{n, l, twoJplus1}, lvIndex - 1});

            // Go the the next line
            getline(elvlc, line);
        }
        elvlc.close();

        // take the appropriate averages (either over j, or over j and l)
        _levelv.clear();
        int n = 1;
        while (n <= _resolvedUpTo)
        {
            for (int l = 0; l < n; l++)
            {
                // average over the j states
                double energy = 0;
                for (int twoJplus1 : twoJplus1range(l))
                    energy += _chiantiLevelv[nljToChiantiIndexm.at({n, l, twoJplus1})].e() * twoJplus1;
                energy /= 4 * l + 2;

                _levelv.emplace_back(n, l, energy);
                _nlToIndexm.insert({{n, l}, _levelv.size() - 1});
            }
            n++;
        }
        while (n <= cNMAX)
        {
            // overage over the l and j states
            double energy = 0;
            for (int l = 0; l < n; l++)
            {
                for (int twoJplus1 : twoJplus1range(l))
                    energy += _chiantiLevelv[nljToChiantiIndexm.at({n, l, twoJplus1})].e() * twoJplus1;
            }
            energy /= 2 * n * n;

            _levelv.emplace_back(n, energy);
            _nlToIndexm.insert({{n, -1}, _levelv.size() - 1});
            n++;
        }
        _numL = _levelv.size();
    }

    void HFromFiles::readTransProbs()
    {
        _chiantiAvv = EMatrix::Zero(_chiantiLevelv.size(), _chiantiLevelv.size());
        ifstream wgfa = IOTools::ifstreamRepoFile(chiantiBaseName + ".wgfa");
        string line;
        getline(wgfa, line);
        while (line.compare(1, 2, "-1"))
        {
            int leftIndex, rightIndex;
            double wavAngstrom, gf, A;
            istringstream(line) >> leftIndex >> rightIndex >> wavAngstrom >> gf >> A;

            // A comment in the cloudy code recommended to do this, as there are apparently some
            // files in the CHIANTI database where the left index represents the upper level of
            // the transition
            int upperIndex = max(leftIndex, rightIndex);
            int lowerIndex = min(leftIndex, rightIndex);

            // Zero means two-photon transition, see CHIANTI user guide.
            if (wavAngstrom > 0) _chiantiAvv(upperIndex - 1, lowerIndex - 1) = A;
            getline(wgfa, line);
        }
        wgfa.close();

        // Total transition rate from each level = sum over each row
        _totalAv = makeAvv().rowwise().sum();
    }

    void HFromFiles::readCollisions()
    {
        // Set up index -> nl translation for the Anderson data. Subtract 1 when using an index
        // from the file!
        for (int n = 0; n <= cNMAX; n++)
            for (int l = 0; l < n; l++) _andersonIndexm1ToNLv.push_back({n, l});

        ifstream h_coll_str = IOTools::ifstreamRepoFile("dat/h_coll_str.dat");
        string line;
        getline(h_coll_str, line);
        getline(h_coll_str, line);

        // For each line, get the initial level index, final level index, and the data
        vector<int> iv, fv;
        vector<Array> qvv;
        while (line.compare(0, 2, "-1"))
        {
            istringstream iss = istringstream(line);
            int u, l;
            Array Upsilonv(8);
            iss >> u >> l;
            for (int t = 0; t < 8; t++) iss >> Upsilonv[t];
            iv.emplace_back(u);
            fv.emplace_back(l);
            qvv.emplace_back(Upsilonv);
            getline(h_coll_str, line);
        }
        h_coll_str.close();

        // These temperatures are not in the file, but they are in the paper, in eV
        Array andersonTemp_eV{{0.5, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0}};
        Array andersonTemp_K = andersonTemp_eV * Constant::EV / Constant::BOLTZMAN;

        // Insert the gathered data into the collision data object
        int numTransitions = qvv.size();
        _qdataAnderson.prepare(andersonTemp_K, numTransitions);
        for (int t = 0; t < numTransitions; t++) _qdataAnderson.insertDataForTransition(qvv[t], iv[t], fv[t]);
        _qdataAnderson.check();
    }

    size_t HFromFiles::index(int n, int l) const
    {
        if (n > _resolvedUpTo)
            return _nlToIndexm.at({n, -1});
        else
            return _nlToIndexm.at({n, l});
    }

    EVector HFromFiles::makeEv() const
    {
        EVector the_ev(_numL);
        for (size_t i = 0; i < _numL; i++) the_ev(i) = _levelv[i].e();
        return the_ev;
    }

    EVector HFromFiles::makeGv() const
    {
        EVector the_gv(_numL);
        for (size_t i = 0; i < _numL; i++) the_gv(i) = _levelv[i].g();
        return the_gv;
    }

    EMatrix HFromFiles::makeAvv() const
    {
        EMatrix the_avv = EMatrix::Zero(_numL, _numL);

        // algorithm: for each transition coefficient available from the data (nlj-resolved),
        // find out the corresponding indices in our Avv matrix, and add the coefficient with
        // the right weight.
        int numLevelsFromData = _chiantiLevelv.size();
        for (int i = 0; i < numLevelsFromData; i++)
        {
            // When combining different initial (upper) levels, the coefficients for decay to a
            // specific final (lower) level need to be averaged using the statistical weights g
            // of the original levels, relative to the collapsed level (note that sum of
            // g_original = g_collapsed). This is based on the assumption that the populations
            // are statistically distributed.
            int iTarget = index(_chiantiLevelv[i].n(), _chiantiLevelv[i].l());
            double weight = _chiantiLevelv[i].g() / _levelv[iTarget].g();

            for (int f = 0; f < numLevelsFromData; f++)
            {
                // When combining different lower levels, the coefficients for decay from a
                // specific upper level simply need to be summed (the total arrival rate into
                // the lower level is the sum of the arrival rates into the other levels). We do
                // not need to check the collapse scheme (fTarget will point to the same level
                // multiple times, if coefficients need to be summed).
                int fTarget = index(_chiantiLevelv[f].n(), _chiantiLevelv[f].l());
                the_avv(iTarget, fTarget) += weight * _chiantiAvv(i, f);
            }
        }
        return the_avv;
    }

    EMatrix HFromFiles::makeExtraAvv() const
    {
        EMatrix the_extra = EMatrix::Zero(_numL, _numL);

        size_t index1s = index(1, 0);

        // Check if the n2 level is resolved, and retrieve the Key, Value pair
        auto index2sIt = _nlToIndexm.find({2, 0});
        // If the level is collapsed, the pair 0,2 won't be found, and the find function will
        // return 'end'. So if the result is not 'end', this means that the resolved 2,0 level was
        bool n2Resolved = index2sIt != _nlToIndexm.end();

        // Hardcode the two-photon decay of 2s to 1s
        double rate = 8.229;
        // If 2nd level is resolved, just put the value in the right place
        if (n2Resolved)
        {
            size_t index2s = index2sIt->second;
            the_extra(index2s, index1s) = rate;
        }
        // If it is collapsed, assume it is well mixed. So we take the average of the two-photon
        // rates from 2s (8.229) and 2p (0).
        else
        {
            size_t index2 = index(2, -1);
            the_extra(index2, index1s) = rate / 4.;  // + 0 * 3 / 4.
        }
        return the_extra;
    }

    array<size_t, 2> HFromFiles::twoPhotonIndices() const
    {
        // If any of the levels is not resolved on l, just return the index of the collapsed level.
        size_t upper = _resolvedUpTo >= 2 ? index(2, 0) : index(2, -1);
        size_t lower = _resolvedUpTo >= 1 ? index(1, 0) : index(1, -1);
        return {upper, lower};
    }

    EMatrix HFromFiles::cvv(const CollisionParameters& cp) const
    {
        double T = cp._t;
        double kT = Constant::BOLTZMAN * T;

        EMatrix the_cvv = EMatrix::Zero(_numL, _numL);

        // Electron contributions (n-changing)
        double ne = cp._sv.ne();
        if (ne > 0)
        {
            // interpolate the collision strength data in temperature space
            Array UpsilonDownv = _qdataAnderson.qv(T);

            // process the interpolated data
            int numTransitions = _qdataAnderson.transitionv().size();
            for (int transitionIndex = 0; transitionIndex < numTransitions; transitionIndex++)
            {
                // get initial and final n,l and the corresponding degeneracy
                int iData = _qdataAnderson.transitionv()[transitionIndex][0];
                int ni = _andersonIndexm1ToNLv[iData - 1][0];
                int li = _andersonIndexm1ToNLv[iData - 1][1];
                double gi = 4 * li + 2;

                int fData = _qdataAnderson.transitionv()[transitionIndex][1];
                int nf = _andersonIndexm1ToNLv[fData - 1][0];
                int lf = _andersonIndexm1ToNLv[fData - 1][1];
                double gf = 4 * li + 2;

                // get initial and final index where the contribution needs to go
                int iTarget = index(ni, li);
                int fTarget = index(nf, lf);

                // upward and downward rates, derived from the interpolated collision strength
                double rateDown = UpsilonDownv[transitionIndex] * 8.6291e-6 / (gi * sqrt(T)) * ne;
                double rateUp = rateDown * gi / gf * exp((_levelv[fTarget].e() - _levelv[iTarget].e()) / kT);

                // analogously to makeAvv(), add these numbers to the transition matrix at the
                // right position and with the right weights
                double iWeight = gi / _levelv[iTarget].g();
                double fWeight = gf / _levelv[fTarget].g();
                the_cvv(iTarget, fTarget) += iWeight * rateDown;
                the_cvv(fTarget, iTarget) += fWeight * rateUp;
            }
        }

        // For the l-resolved levels, get l-changing collision rates (through proton collisions)
        double np = cp._sv.np();
        if (np > 0)
        {
            for (int n = 0; n <= _resolvedUpTo; n++)
            {
                // The formulae of PS64 are implemented separately. The result is indexed on li,
                // lf.
                EMatrix qvv = PS64CollisionRateCoeff(n, T, np);
                // Fill in the collision rates for all combinations of li lf
                for (int li = 0; li < n; li++)
                {
                    size_t i = index(n, li);
                    for (int lf = 0; lf < n; lf++)
                    {
                        size_t f = index(n, lf);
                        // None of the previous contributions should have been l-changing
                        assert(the_cvv(i, f) == 0);
                        the_cvv(i, f) += qvv(li, lf) * np;
                    }
                }
            }
        }
        // Make all negative entries 0
        return the_cvv.array().max(0);
    }

    EMatrix HFromFiles::PS64CollisionRateCoeff(int n, double T, double np) const
    {
        EMatrix q_li_lf_goingUp = EMatrix::Zero(n, n);
        EMatrix q_li_lf_goingDown = EMatrix::Zero(n, n);

        // We will apply PS64 eq 43. Keeping eq 38 in mind, we can find the partial rates one by
        // one, applying detailed balance at each step. Note however that the results are
        // different depending on which side (l = 0 or l = n - 1) you start. We calculate the
        // coefficients in both ways and then take the average.

        // mu is the reduced mass of the system of the colliding particles
        constexpr double muOverm =
            Constant::PROTONMASS * Constant::HMASS / (Constant::PROTONMASS + Constant::HMASS) / Constant::ELECTRONMASS;

        const double qnlFactor = 9.93e-6 * sqrt(muOverm / T);
        const int n2 = n * n;

        auto ps64eq43 = [&](int l) -> double {
            // eq 44: Z is charge of the colliding particle, z that of the nucleus
            double D_nl = 6 * n2 * (n2 - l * l - l - 1);

            // eq 45,46: take the smallest of the two, since Rc represents a cutoff value that
            // prevented divergence in the calculations of PS64
            size_t i = index(n, l);
            double tau2 = 1. / _totalAv(i) / _totalAv(i);
            double twoLog10Rc = min(10.95 + log10(T * tau2 / muOverm), 1.68 + log10(T / np));

            // eq 43
            double q_nl = qnlFactor * D_nl * (11.54 + log10(T / D_nl / muOverm) + twoLog10Rc);
            return max(0., q_nl);
        };

        double qUp = 0;
        double qDown = 0;
        // The value for q(l = n - 1, l = n - 2) is filled in in the last iteration
        for (int l = 0; l < n - 1; l++)
        {
            double qBoth = ps64eq43(l);
            qUp = qBoth - qDown;
            q_li_lf_goingUp(l, l + 1) = qUp;

            // will also be used in loop for l+1:
            qDown = qUp * (2. * l + 1.) / (2. * l + 3.);
            q_li_lf_goingUp(l + 1, l) = qDown;
        }
        qUp = 0;
        qDown = 0;
        for (int l = n - 1; l > 0; l--)
        {
            double qBoth = ps64eq43(l);
            qDown = qBoth - qUp;
            q_li_lf_goingDown(l, l - 1) = qDown;

            // will also be used in loop for l-1:
            qUp = qDown * (2. * l + 1.) / (2. * l - 1.);
            q_li_lf_goingDown(l - 1, l) = qUp;
        }
        const EMatrix result = (q_li_lf_goingUp + q_li_lf_goingDown) / 2.;
        return result.array().max(0);
    }
}
