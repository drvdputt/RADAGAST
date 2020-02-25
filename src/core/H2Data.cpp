#include "H2Data.hpp"
#include "CollisionParameters.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"
#include <regex>

using namespace std;

H2Data::H2Data(int maxJ, int maxV) : LevelCoefficients(2 * Constant::HMASS), _maxJ{maxJ}, _maxV{maxV}
{
    readLevels();
    readTransProbs();

    // Set these data members of the LevelCoefficients parent class
    setConstants(makeEv(), makeGv(), _avv, EMatrix::Zero(_numL, _numL));

    readCollisions();
    readDissProbs();
    readDirectDissociation();

    // Now that all the leves are in order, precalculate the g-bar rates
    precalcGBarKvv();

    if (Options::h2data_plotLevelMatrices)
    {
        ofstream avvOut = IOTools::ofstreamFile("h2/einsteinA.dat");
        avvOut << _avv << '\n';
        avvOut.close();

        ofstream lvlOut = IOTools::ofstreamFile("h2/levels.dat");
        lvlOut << "# eState\tv\tJ\tE\n";
        for (const auto& l : _levelv)
            lvlOut << static_cast<int>(l.eState()) << "\t" << l.v() << "\t" << l.j() << "\t" << l.e() << '\n';
        lvlOut.close();
    }
}

size_t H2Data::indexOutput(ElectronicState eState, int j, int v) const
{
    // This cast is safe, as the default underlying type of an enum is int.
    return _ejvToIndexm.at({static_cast<int>(eState), j, v});
}

int H2Data::indexFind(ElectronicState eState, int j, int v) const
{
    auto iter = _ejvToIndexm.find({static_cast<int>(eState), j, v});
    if (iter == _ejvToIndexm.end())
        return -1;
    else
        return (*iter).second;
}

void H2Data::readLevels()
{
    // Expand levelv with the levels listed in these files
    readLevelFile("dat/h2/energy_X.dat", ElectronicState::X);
    _startOfExcitedIndices = _levelv.size();
    if (Options::h2data_numExcitedLevels >= 1)
    {
        // Lyman
        readLevelFile("dat/h2/energy_B.dat", ElectronicState::B);
        _bB = true;
    }
    if (Options::h2data_numExcitedLevels >= 2)
    {
        // Werner
        readLevelFile("dat/h2/energy_C_plus.dat", ElectronicState::Cplus);
        _bCplus = true;
    }
    if (Options::h2data_numExcitedLevels >= 3)
    {
        readLevelFile("dat/h2/energy_C_minus.dat", ElectronicState::Cminus);
        _bCminus = true;
    }
    _numL = _levelv.size();
}

void H2Data::readTransProbs()
{
    _avv = EMatrix::Zero(_numL, _numL);

    // Radiative transitions within ground stated. Fills _avv with data from Wolniewicz
    // (1998)
    readTransProbFile("dat/h2/transprob_X.dat", ElectronicState::X, ElectronicState::X);
    if (_bB) readTransProbFile("dat/h2/transprob_B.dat", ElectronicState::B, ElectronicState::X);
    if (_bCplus) readTransProbFile("dat/h2/transprob_C_plus.dat", ElectronicState::Cplus, ElectronicState::X);
    if (_bCminus) readTransProbFile("dat/h2/transprob_C_minus.dat", ElectronicState::Cminus, ElectronicState::X);
}

void H2Data::readDissProbs()
{
    _dissProbv = EVector::Zero(_numL);
    _dissKinEv = EVector::Zero(_numL);

    if (_bB) readDissProbFile("dat/h2/dissprob_B.dat", ElectronicState::B);
    if (_bCplus) readDissProbFile("dat/h2/dissprob_C_plus.dat", ElectronicState::Cplus);
    if (_bCminus) readDissProbFile("dat/h2/dissprob_C_minus.dat", ElectronicState::Cminus);
}

void H2Data::readCollisions()
{
    _qdataPerPartner.resize(_numPartners);
    _hasQdata.resize(_numPartners, EMatrix_bool::Zero(_numL, _numL));

    // Collisions with H, from Lique 2015
    readCollisionFile("dat/h2/coll_rates_H_15.dat", H0);

    // Collision rates with H2, from Lee 2008
    readCollisionFile("dat/h2/coll_rates_H2ortho_ORNL.dat", H2ORTHO);
    readCollisionFile("dat/h2/coll_rates_H2para_ORNL.dat", H2PARA);

    // Collision rates with H+, from Gerlich 1990
    readCollisionFile("dat/h2/coll_rates_Hp.dat", HPLUS);
}

void H2Data::readDirectDissociation()
{
    _dissociationCrossSectionv.resize(_levelv.size());
    ifstream cont_diss = IOTools::ifstreamRepoFile("dat/h2/cont_diss.dat");
    string line;

    // Main file traversal loop
    while (getline(cont_diss, line))
    {
        // Find the start of a cross section data block
        if (line.at(0) == '#' && line.at(1) == '!')
        {
            // Get the three lines starting with "#!". I don't think I need the
            // 2nd and 3rd lines though.
            string line1, line2;
            getline(cont_diss, line1);
            getline(cont_diss, line2);

            // Parse the first line to get the quantum numbers
            int nei, nef, vi, ji;
            string cleanline = regex_replace(line, regex("[^0-9.]+"), " ");
            istringstream(cleanline) >> nei >> nef >> vi >> ji;

            // Skip to the next #! block if we are not treating this level (There
            // should be no data for states other than X, but let's ignore nei
            // anyway. .
            if (nei > 0 || !validJV(ji, vi)) continue;  // This will never be used;

            // Read the data
            string dataline;
            vector<double> frequencyv, crossSectionv;

            while (getline(cont_diss, dataline))
            {
                // Parse the data line:	Split by comma
                istringstream iss(dataline);
                string s0, s1;
                getline(iss, s0, ',');
                getline(iss, s1);

                // String to double
                double energy_invcm = stod(s0);
                double crossSection_ang2 = stod(s1);

                // Convert to the correct units and add
                frequencyv.emplace_back(Constant::LIGHT * energy_invcm);
                crossSectionv.emplace_back(crossSection_ang2 * Constant::ANG_CM * Constant::ANG_CM);

                if (cont_diss.peek() == '#') break;
            }
            // When this loop exits, we should be at the start of the next comment
            // block, or at the end of the file.

            // Add the new cross section at the correct level index.
            int index = indexFind(ElectronicState::X, ji, vi);
            if (index > -1)
                _dissociationCrossSectionv[index].emplace_back(Array(frequencyv.data(), frequencyv.size()),
                                                               Array(crossSectionv.data(), crossSectionv.size()));
        }
    }
    // Remember which levels have cross sections, for easy iterating
    _levelsWithCrossSectionv.reserve(_numL);
    for (size_t i = 0; i < _numL; i++)
        if (!_dissociationCrossSectionv[i].empty()) _levelsWithCrossSectionv.emplace_back(i);
    _levelsWithCrossSectionv.shrink_to_fit();
}

void H2Data::readLevelFile(const string& repoFile, ElectronicState eState)
{
    ifstream energy = IOTools::ifstreamRepoFile(repoFile);

    string line;
    getline(energy, line);
    int y, m, d;
    istringstream(line) >> y >> m >> d;
    Error::equalCheck("magic number y", y, 2);
    Error::equalCheck("magic number m", m, 4);
    Error::equalCheck("magic number d", d, 29);

    // Start reading in the V, J and E (cm-1) of the levels
    vector<int> vv, jv;
    vector<double> ev;
    int counter = 0;
    while (getline(energy, line))
    {
        auto iss = istringstream(line);
        if (line.front() == '#')
        {
            string word1, word2;
            iss >> word1 >> word2;
            continue;  // Do not read the 'extra' levels to better match with cloudy
        }
        int v, j;
        iss >> v >> j;
        double k;  // cm-1
        iss >> k;

        // Add the level if within the requested J,v limits. Convert 1/wavelength [cm-1]
        // to ergs by using E = hc / lambda.
        if (validJV(j, v))
        {
            counter++;
            addLevel(eState, j, v, Constant::PLANCKLIGHT * k);
        }
    }
    DEBUG("Read in " << counter << " levels from " << repoFile << '\n');
}

void H2Data::readTransProbFile(const string& repoFile, ElectronicState upperE, ElectronicState lowerE)
{
    ifstream transprob = IOTools::ifstreamRepoFile(repoFile);

    string line;
    getline(transprob, line);
    int y, m, d;
    istringstream(line) >> y >> m >> d;
    Error::equalCheck("magic number y", y, 2);
    Error::equalCheck("magic number m", m, 4);
    Error::equalCheck("magic number d", d, 29);

    size_t counter = 0;
    while (getline(transprob, line))
    {
        if (line.empty() || line.at(0) == '#') continue;

        int EU, VU, JU, EL, VL, JL;
        double A;  // Unit s-1
        istringstream(line) >> EU >> VU >> JU >> EL >> VL >> JL >> A;

        // Check if the entry in the file and the expected Estate-changing transition
        // correspond (just a double check on the enum scheme)
        Error::equalCheck("EU and cast of upperE", EU, static_cast<int>(upperE));
        Error::equalCheck("EL and cast of lowerE", EL, static_cast<int>(lowerE));

        // If within the J,v limits, fill in the coefficient.
        if (validJV(JU, VU) && validJV(JL, VL))
        {
            int upperIndex = indexFind(upperE, JU, VU);
            int lowerIndex = indexFind(lowerE, JL, VL);
            if (upperIndex > -1 && lowerIndex > -1)
            {
                _avv(upperIndex, lowerIndex) = A;
                counter++;
            }
        }
    }
    DEBUG("Read in " << counter << " Einstein A coefficients from " << repoFile << '\n');
}

void H2Data::readDissProbFile(const string& repoFile, ElectronicState eState)
{
    ifstream dissprobs = IOTools::ifstreamRepoFile(repoFile);
    int y, m, d;
    IOTools::istringstreamNextLine(dissprobs) >> y >> m >> d;
    Error::equalCheck("magic y", y, 3);
    Error::equalCheck("magic m", m, 2);
    Error::equalCheck("magic d", d, 11);

    string line;
    int counter = 0;
    while (getline(dissprobs, line))
    {
        if (line.empty() || line.at(0) == '#') continue;

        int v, j;
        // diss prob in s-1, kin energy in eV
        double diss, kin;
        istringstream(line) >> v >> j >> diss >> kin;

        if (!validJV(j, v)) continue;

        int index = indexFind(eState, j, v);
        if (index > -1)
        {
            _dissProbv[index] = diss;
            _dissKinEv[index] = kin / Constant::ERG_EV;
            counter++;
        }
    }
    DEBUG("Read in " << counter << " dissociation rates from " << repoFile << '\n');
}

void H2Data::readCollisionFile(const string& repoFile, CollisionPartner iPartner)
{
    ifstream coll_rates = IOTools::ifstreamRepoFile(repoFile);

    int m;
    IOTools::istringstreamNextLine(coll_rates) >> m;
    Error::equalCheck(repoFile + "magic number", m, 110416);

    // Read until the temperature line
    string line;
    while (getline(coll_rates, line))
        if (line.at(0) != '#') break;

    // The number of dots gives us the number of columns
    size_t numTemperatures = count(begin(line), end(line), '.');

    // Read in the temperatures
    Array temperaturev(numTemperatures);
    istringstream issTemperatures(line);
    for (size_t i = 0; i < numTemperatures; i++) issTemperatures >> temperaturev[i];

    // For each line, get the initial level index, final level index, and the data
    vector<int> iv, fv;
    vector<Array> qvv;
    while (getline(coll_rates, line))
    {
        if (line.empty() || line.at(0) == '#') continue;

        istringstream issCollRates(line);
        int VU, JU, VL, JL;
        issCollRates >> VU >> JU >> VL >> JL;

        // If we are not treating this level (the J or the v is out of range), skip.
        if (!validJV(JU, VU) || !validJV(JL, VL)) continue;

        int upperIndex = indexFind(ElectronicState::X, JU, VU);
        int lowerIndex = indexFind(ElectronicState::X, JL, VL);
        if (upperIndex > -1 && lowerIndex > -1)
        {
            iv.emplace_back(upperIndex);
            fv.emplace_back(lowerIndex);
            Array qForEachTv(numTemperatures);
            for (size_t i = 0; i < numTemperatures; i++) issCollRates >> qForEachTv[i];
            qvv.emplace_back(qForEachTv);
        }
        // else the data is unused, and we go to the next line
    }

    size_t numTransitions = qvv.size();
    _qdataPerPartner[iPartner].prepare(temperaturev, numTransitions);
    for (size_t t = 0; t < numTransitions; t++)
    {
        _qdataPerPartner[iPartner].insertDataForTransition(qvv[t], iv[t], fv[t]);
        _hasQdata[iPartner](iv[t], fv[t]) = true;
    }
    _qdataPerPartner[iPartner].check();
    DEBUG("Read in " << _qdataPerPartner[iPartner].transitionv().size() << " collision coefficients from " << repoFile
                     << '\n');
}

EVector H2Data::makeEv() const
{
    EVector the_ev(_numL);
    for (size_t i = 0; i < _numL; i++) the_ev(i) = _levelv[i].e();
    return the_ev;
}

EVector H2Data::makeGv() const
{
    EVector the_gv(_numL);
    for (size_t i = 0; i < _numL; i++) the_gv(i) = _levelv[i].g();
    return the_gv;
}

EMatrix H2Data::cvv(const CollisionParameters& cp) const
{
    double T = cp._t;
    EMatrix the_cvv{EMatrix::Zero(_numL, _numL)};

    // H-H2 collisions
    double nH = cp._sv.nH();
    addToCvv(the_cvv, T, H0, nH);

    // H2-H2 collisions.
    double nH2 = cp._sv.nH2();
    addToCvv(the_cvv, T, H2ORTHO, cp._orthoH2 * nH2);
    addToCvv(the_cvv, T, H2PARA, (1. - cp._orthoH2) * nH2);

    // TODO: test if these are loaded correctly (compare with graph in Lee paper)

    // H+-H2 collisions
    double np = cp._sv.np();
    addToCvv(the_cvv, T, HPLUS, np);
    return the_cvv;
}

double H2Data::directDissociationCrossSection(double nu, int j, int v) const
{
    return directDissociationCrossSection(nu, indexOutput(ElectronicState::X, j, v));
}

double H2Data::directDissociationCrossSection(double nu, size_t index) const
{
    double sigma{0.};
    // Evaluate all the cross sections for this level at this frequency
    for (const Spectrum& cs : _dissociationCrossSectionv[index]) sigma += cs.evaluate(nu);
    return sigma;
}

const vector<Spectrum>& H2Data::directDissociationCrossSections(size_t index) const
{
    return _dissociationCrossSectionv[index];
}

bool H2Data::validJV(int J, int v) const
{
    return J <= _maxJ && v <= _maxV;
}

void H2Data::addToCvv(EMatrix& the_cvv, double T, CollisionPartner iPartner, double nPartner) const
{
    if (nPartner <= 0) return;

    const CollisionData& qdata = _qdataPerPartner[iPartner];

    /* Find the grid point to the right of (>=) the requested log-temperature. (Returns last
	   point if T > Tmax). We will naively extrapolate for points outside the range, and cut
	   off the result at 0 if it becomes negative. */
    const Array& temperaturev = qdata.temperaturev();
    int iRight = TemplatedUtils::index(T, temperaturev);
    iRight = max(iRight, 1);
    int iLeft = iRight - 1;
    double tRight = temperaturev[iRight];
    double tLeft = temperaturev[iLeft];

    double kT = T * Constant::BOLTZMAN;
    const vector<array<int, 2>>& transitionv = qdata.transitionv();
    for (int transitionIndex = 0; transitionIndex < transitionv.size(); transitionIndex++)
    {
        // Interpolate data points naively for temperatures left and right of T
        double qLeft = qdata.q(iLeft, transitionIndex);
        double qRight = qdata.q(iRight, transitionIndex);
        double q = TemplatedUtils::interpolateLinear(T, tLeft, tRight, qLeft, qRight);
        // Make zero if interpolation makes this negative
        q = max(0., q);

        // The levels involved in this transition
        int i = transitionv[transitionIndex][0];
        int f = transitionv[transitionIndex][1];
        double Cif = q * nPartner;
        the_cvv(i, f) += Cif;
        the_cvv(f, i) += otherDirectionC(Cif, i, f, kT);
    }
    addGBarCvv(the_cvv, kT, iPartner, nPartner);
}

void H2Data::precalcGBarKvv()
{
    _gBarKvvPerPartner.resize(_numPartners, EMatrix::Zero(_startOfExcitedIndices, _startOfExcitedIndices));
    for (int p = 0; p < _numPartners; p++)
    {
        for (int i = 0; i < _startOfExcitedIndices; i++)
        {
            for (int j = i + 1; j < _startOfExcitedIndices; j++)
            {
                // apply eq 1 from Shaw et al. 2005. sigma is transition energy
                // in cm-1
                double sigma = (_levelv[i].e() - _levelv[j].e()) / Constant::PLANCKLIGHT;
                double y0 = _gbarcoll[p][0];
                double a = _gbarcoll[p][1];
                double b = _gbarcoll[p][2];
                double logk = y0 + a * pow(max(abs(sigma), 100.), b);
                double kDown = pow(10., logk);

                // Only fill in downward rates. Upward rates will depend on
                // temperature using otherDirectionC.
                if (sigma > 0)  // e_i > e_j
                    _gBarKvvPerPartner[p](i, j) = kDown;
                else if (sigma < 0)  // e_j > e_i
                    _gBarKvvPerPartner[p](j, i) = kDown;
            }
        }
    }
}

void H2Data::addGBarCvv(EMatrix& the_cvv, double kT, CollisionPartner iPartner, double nPartner) const
{
    for (int i = 0; i < _startOfExcitedIndices; i++)
    {
        for (int j = i + 1; j < _startOfExcitedIndices; j++)
        {
            // If we already have data for this collision partner and level, skip
            if (_hasQdata[iPartner](i, j)) continue;

            // Determine which is upper and which is lower (I don't know of j > i
            // guarantees that e_j > e_i)
            int u, l;
            if (_levelv[i].e() > _levelv[j].e())
            {
                u = i;
                l = j;
            }
            else
            {
                u = j;
                l = i;
            }
            double Cul = nPartner * _gBarKvvPerPartner[iPartner](u, l);
            the_cvv(u, l) += Cul;
            the_cvv(l, u) += otherDirectionC(Cul, u, l, kT);
        }
    }
}

double H2Data::otherDirectionC(double Cif, int i, int f, double kT) const
{
    double Cfi =
        Cif * _levelv[i].g() / static_cast<double>(_levelv[f].g()) * exp((_levelv[f].e() - _levelv[i].e()) / kT);
    return Cfi;
}

EVector H2Data::formationDistribution() const
{
    double kTf = Constant::BOLTZMAN * 5.e4;

    EVector fv = EVector(_startOfExcitedIndices);
    for (int i = 0; i < _startOfExcitedIndices; i++)
    {
        int gj = (_levelv[i].j() % 2) ? 3 : 1;
        fv[i] = gj * (1 + _levelv[i].v()) * exp(-_levelv[i].e() / kTf);
    }
    return fv / fv.sum();
}
