#include "GrainPhotoelectricCalculator.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "RadiationFieldTools.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "WeingartnerDraine2001.hpp"
#include <cmath>
#include <iomanip>

using namespace std;

GrainPhotoelectricCalculator::GrainPhotoelectricCalculator(const Array& sizev, double workFunction, bool carOrSil)
    : _workFunction{workFunction}, _carOrSil{carOrSil}, _sizev{sizev}
{
    _y1Cache.resize(_sizev.size());
    _eStickPositiveCache.resize(_sizev.size());
    _eStickNegativeCache.resize(_sizev.size());
    for (size_t i = 0; i < _sizev.size(); i++)
    {
        _y1Cache[i] = WD01::y1(_sizev[i]);
        _eStickPositiveCache[i] = WD01::estick_positive(_sizev[i]);
        _eStickNegativeCache[i] = WD01::estick_negative(_sizev[i]);
    }
}

GrainPhotoelectricCalculator::Locals::Locals(const Spectrum* specificIntensity, double T, double ne,
                                             const std::vector<int>& chargev, const Array& densityv, const Array& massv)
    : _specificIntensity(specificIntensity), _T(T), _ne(ne), _chargev(chargev), _densityv(densityv), _massv(massv),
      _integrationWorkspace(specificIntensity->numPoints())
{}

int GrainPhotoelectricCalculator::minimumCharge(int i) const
{
    double Uait = autoIonizationThreshold(i);
    return WD01::minimumCharge(_sizev[i], Uait);
}

void GrainPhotoelectricCalculator::calculateChargeDistribution(int i, Locals& env, const Array& Qabsv,
                                                               ChargeDistribution& cd) const
{
    double a = _sizev[i];

    // Express a in angstroms
    double aA = a / Constant::ANG_CM;

    // Highest possible energy of a photon
    const Array& frequencyv = env._specificIntensity->frequencyv();
    double hnumax = Constant::PLANCK * frequencyv[frequencyv.size() - 1];

    // The maximum charge is one more than the highest charge which still allows ionization by
    // photons of hnumax.
    int resultZmax = floor(((hnumax - _workFunction) * Constant::ERG_EV / 14.4 * aA + .5 - .3 / aA) / (1 + .3 / aA));

    // The minimum charge is the most negative charge for which autoionization does not occur
    int resultZmin = minimumCharge(i);

    // The few cases I've seen this happen, Zmin is always 0 and Zmax is always -1. Since Zmax only
    // depends on the maximum photon energy, it should be fine to increase it up to Zmin in these
    // edge cases.
    if (resultZmax < resultZmin) resultZmax = resultZmin;

    // These two functions determine the up and down rates for the detailed balance

    // The rate at which the grain moves out of charge Z, in the positive direction. Contributions
    // by photoelectric ejection of electron, and collisions with positive particles. Collisions
    // with multiply charged particles are treated as if they can only change the charge in steps
    // of 1. This is technically incorrect, but the charging rate due to positive ions is very
    // small anyway.
    auto chargeUpRate = [&](int z) -> double {
        double Jtotal{0};
        Jtotal += emissionRate(i, z, env, Qabsv);
        for (size_t j = 0; j < env._chargev.size(); j++)
        {
            if (env._chargev[j] > 0)
                Jtotal += collisionalChargingRate(i, env._T, z, env._chargev[j], env._massv[j], env._densityv[j]);
        }
        return Jtotal;
    };
    // The rate at which the grain moves out of charge Z, in the negative direction. Sums
    // collisions over negative particles.
    auto chargeDownRate = [&](int z) -> double {
        double Jtotal{0};
        for (size_t j = 0; j < env._chargev.size(); j++)
        {
            if (env._chargev[j] < 0)
                Jtotal += collisionalChargingRate(i, env._T, z, env._chargev[j], env._massv[j], env._densityv[j]);
        }
        return Jtotal;
    };
    // Limit the number of charges here
    cd.calculateDetailedBalance(chargeUpRate, chargeDownRate, resultZmin, resultZmax,
                                Options::grainphotoelectriccalculator_maxcharges);
    DEBUG("grain" << i << " average charge " << cd.average() << '\n');
}

void GrainPhotoelectricCalculator::getPET_PDT_Emin(int i, int Z, double& pet, double& pdt, double& Emin) const
{
    Emin = WD01::eMin(_sizev[i], Z);
    double ip_v = ionizationPotential(i, Z);
    // WD01 eq 6
    pet = Z >= -1 ? ip_v : ip_v + Emin;
    // Photodetachment, WD01 eq 18 (using EA(Z + 1, a) = IP(Z, a) if Z < 0)
    pdt = ip_v + Emin;
}

double
GrainPhotoelectricCalculator::photoelectricIntegrationLoop(Locals& env, const Array& Qabsv, double pet,
                                                           const std::function<double(double hnuDiff)>* f_hnuDiff) const
{
    const Array& frequencyv = env._specificIntensity->frequencyv();
    const Array& specificIntensityv = env._specificIntensity->valuev();

    // yield is zero by definition below photoelectric threshold
    double nu_pet = pet / Constant::PLANCK;
    int i = frequencyv.size() - 1;
    while (frequencyv[i] > nu_pet && i > 0)
    {
        double hnu = Constant::PLANCK * frequencyv[i];
        env._integrationWorkspace[i] = Qabsv[i] * specificIntensityv[i] / hnu;
        if (f_hnuDiff) env._integrationWorkspace[i] *= (*f_hnuDiff)(hnu - pet);
        i--;
    }
    // <function unit> sr-1 cm-2 s-1
    double integral =
        TemplatedUtils::integrate<double>(frequencyv, env._integrationWorkspace, max(0, i), frequencyv.size() - 1);
    // <function unit> cm-2 s-1
    return Constant::FPI * integral;
}

double GrainPhotoelectricCalculator::photodetachmentIntegrationLoop(int Z, Locals& env, double pdt,
                                                                    const double* calcEnergyWithThisEmin) const
{
    const Array& frequencyv = env._specificIntensity->frequencyv();
    const Array& specificIntensityv = env._specificIntensity->valuev();

    // no effect below photodetachment threshold
    double nu_pdt = pdt / Constant::PLANCK;
    int i = frequencyv.size() - 1;
    while (frequencyv[i] > nu_pdt && i > 0)
    {
        double hnu = Constant::PLANCK * frequencyv[i];
        double hnuDiff = hnu - pdt;
        // <function unit> / time / angle
        env._integrationWorkspace[i] = WD01::sigmaPDT(Z, hnuDiff) * specificIntensityv[i] / hnu;
        if (calcEnergyWithThisEmin) env._integrationWorkspace[i] *= hnuDiff + *calcEnergyWithThisEmin;
        i--;
    }
    // <erg optional> s-1 sr-1
    double integral =
        TemplatedUtils::integrate<double>(frequencyv, env._integrationWorkspace, max(0, i), frequencyv.size() - 1);
    // <erg optional> s-1
    return Constant::FPI * integral;
}

double GrainPhotoelectricCalculator::heatingRateAZ(int i, int Z, Locals& env, const Array& Qabsv) const
{
    double a = _sizev[i];
    const double e2_a = Constant::ESQUARE / a;

    // It's cheaper to calculate these together (and safer, less code duplication), so I made this
    // funtion, and pass these values as arguments when needed.
    double pet, pdt, Emin;
    getPET_PDT_Emin(i, Z, pet, pdt, Emin);

    // WD01 text between eq 10 and 11
    double Elow = Z < 0 ? Emin : -(Z + 1) * e2_a;

    std::function<double(double)> yieldTimesAverageEnergy = [&](double hnuDiff) -> double {
        double Y = photoelectricYield(i, Z, hnuDiff, Emin);

        double Ehigh = Z < 0 ? hnuDiff + Emin : hnuDiff;
        // The integral over the electron energy distribution (integral E f(E) dE), over the energy
        // range for which electrons can escape
        double IntE = WD01::energyIntegral(Elow, Ehigh, Emin);
        // Divide by (integral f(E) dE) over the same range, which normalizes the above value -->
        // IntE / y2 gives an average energy
        double y2 = WD01::escapingFraction(Z, Elow, Ehigh);

        return Y * IntE / y2;
    };
    double heatingRatePE =
        Constant::PI * a * a * photoelectricIntegrationLoop(env, Qabsv, pet, &yieldTimesAverageEnergy);

    double heatingRatePD = 0;
    if (Z < 0) heatingRatePD = photodetachmentIntegrationLoop(Z, env, pdt, &Emin);
    return heatingRatePE + heatingRatePD;
}

double GrainPhotoelectricCalculator::heatingRateA(int i, Locals& env, const Array& Qabsv,
                                                  const ChargeDistribution& cd) const
{
    double totalHeatingForGrainSize = 0;
    for (int Z = cd.zmin(); Z <= cd.zmax(); Z++)
    {
        double fZz = cd.value(Z);
        if (!isfinite(fZz)) Error::runtime("nan in charge distribution");
        double heatAZ = heatingRateAZ(i, Z, env, Qabsv);

        // Fraction of grains in this charge state * heating by a single particle of charge Z.
        totalHeatingForGrainSize += fZz * heatAZ;
    }
    return totalHeatingForGrainSize;
}

double GrainPhotoelectricCalculator::emissionRate(int i, int Z, Locals& env, const Array& Qabsv) const
{
    double a = _sizev[i];
    double pet, pdt, Emin;
    getPET_PDT_Emin(i, Z, pet, pdt, Emin);

    function<double(double)> yieldf = [&](double hnuDiff) { return photoelectricYield(i, Z, hnuDiff, Emin); };
    double emissionRatePE = Constant::PI * a * a * photoelectricIntegrationLoop(env, Qabsv, pet, &yieldf);
    double emissionRatePD = 0;
    if (Z < 0) emissionRatePD = photodetachmentIntegrationLoop(Z, env, pdt, nullptr);

    return emissionRatePE + emissionRatePD;
}

double GrainPhotoelectricCalculator::collisionalChargingRate(int i, double gasT, int Z, int particleCharge,
                                                             double particleMass, double particleDensity) const
{
    double a = _sizev[i];
    double kT = Constant::BOLTZMAN * gasT;

    // WD eq 26: akT / q^2 = akT / e^2 / z^2
    double tau = a * kT / Constant::ESQUARE / particleCharge / particleCharge;
    // Ze / q = Z / z
    double ksi = static_cast<double>(Z) / static_cast<double>(particleCharge);

    double Jtilde;
    if (ksi < 0)
    {
        Jtilde = (1. - ksi / tau) * (1. + sqrt(2. / (tau - 2. * ksi)));
    }
    else if (ksi > 0)
    {
        double toSquare = 1. + 1. / sqrt(4. * tau + 3. * ksi);
        Jtilde = toSquare * toSquare * exp(-WD01::thetaKsi(ksi) / tau);
    }
    else
    {
        Jtilde = 1. + sqrt(Constant::PI / 2. / tau);
    }
    return particleDensity * stickingCoefficient(i, Z, particleCharge) * sqrt(8. * kT * Constant::PI / particleMass) * a
           * a * Jtilde;
}

double GrainPhotoelectricCalculator::recombinationCoolingRate(int i, const Locals& env,
                                                              const ChargeDistribution& cd) const
{
    double a = _sizev[i];

    // Calculates WD01 equation 42
    double kT = Constant::BOLTZMAN * env._T;
    double eightkT3DivPi = 8 * kT * kT * kT / Constant::PI;

    // For every charged collision partner (effectively electrons and protons), add the
    // contributions for each possible grain charge.
    double particleSum = 0;
    for (size_t j = 0; j < env._chargev.size(); j++)
    {
        int z_j = env._chargev[j];
        if (z_j)
        {
            // tau = akT / q^2 (WD01 eq 26)
            double tau = a * kT / z_j / z_j / Constant::ESQUARE;
            double Zsum = cd.sumOverCharge([&](int zGrain) {
                double ksi = zGrain / static_cast<double>(z_j);
                return stickingCoefficient(i, zGrain, z_j) * WD01::lambdaTilde(tau, ksi);
            });
            particleSum += env._densityv[j] * sqrt(eightkT3DivPi / env._massv[j]) * Zsum;
        }
    }

    // The second term of equation 42: autoionization of grains with the most negative charge
    // inhibits the cooling of the gas. EA(Zmin) = IP(Zmin-1) because IP(Z) = EA(Z+1)
    double secondTerm = 0;
    // This term is only included when the population of the maximally negative grain charge
    // minimumCharge is significant. If it is not siginicant, then fZ will not cover minimumCharge,
    // (and Zmin > minimumCharge).
    int zmin = cd.zmin();
    if (zmin == minimumCharge(i))
        secondTerm = cd.value(zmin) * collisionalChargingRate(i, env._T, zmin, -1, Constant::ELECTRONMASS, env._ne)
                     * ionizationPotential(i, zmin - 1);

    return Constant::PI * a * a * particleSum * kT + secondTerm;
}

double GrainPhotoelectricCalculator::gasGrainCollisionCooling(int i, const Locals& env, const ChargeDistribution& cd,
                                                              double Tgrain, bool forGrain) const
{
    double a = _sizev[i];
    double kT = env._T * Constant::BOLTZMAN;
    double kTgrain = Tgrain * Constant::BOLTZMAN;

    double lambdaG = cd.sumOverCharge([&](int zGrain) {
        double lambdaG_for_this_z = 0;
        // Not entirely sure if this is the right potential
        double Ug = ionizationPotential(i, zGrain);
        double Vg = sqrt(Constant::ESQUARE) * Ug;
        for (size_t j = 0; j < env._massv.size(); j++)
        {
            // Dimensionless
            int z_j = env._chargev[j];
            double ZVg = z_j * Vg;
            double psi = ZVg / kT;
            double eta = psi <= 0 ? 1 - psi : exp(-psi);
            double ksi = psi <= 0 ? 1 - psi / 2 : (1 + psi / 2) * exp(-psi);
            // incoming collision
            double total = 2 * kT * ksi;

            // energy carried away by evaporation of the neutralized/thermalized version of the
            // colliding particle (does not count for electrons, see text below eq 29)
            if (env._chargev[j] >= 0) total -= 2 * kTgrain * eta;

            if (forGrain)
            {
                // A charged particle will slow down or speed up before it hits a charged grain
                total -= ZVg * eta;

                // when protons charge the grain, this means they recombine on the surface (not
                // included in chemical network yet though). The recombination energy needs to go
                // somewhere, so the grain gets extra heat. TODO: don't hardcode these things like
                // this. For now this is not a problem because the proton is the only particle with
                // charge 1. I'm pretty sure using 1 rydberg (= 13.6 eV = ionization potential) is
                // fine here.
                if (env._chargev[j] == 1) total += Constant::RYDBERG * eta;
            }

            double S = 0.;
            if (z_j)
            {
                S = stickingCoefficient(i, zGrain, env._chargev[j]);
            }
            else
            {
                // TODO: for neutral particles, S needs to be the accomodation coefficient instead
                // of the sticking coefficient. Baldwin states that it should be mM / (m^2 + M^2)
                // and references Draine (1978), with M the mass of a typical atom in the grain.
                double m = env._massv[j];
                double M = WD01::atomMass(_carOrSil);
                S = m * M / (m * m + M * M);
            }

            // cm s-1
            double vbar = sqrt(8 * kT / Constant::PI / env._massv[j]);

            // cm-3 * cm s-1 * erg = cm-2 s-1 erg
            lambdaG_for_this_z += env._densityv[j] * vbar * S * total;
        }
        return lambdaG_for_this_z;
    });
    // Remember that 1991-Baldwin gives the rate 'per unit projected area'
    return Constant::PI * a * a * lambdaG;
}

double GrainPhotoelectricCalculator::yieldFunctionTest() const
{
    // Parameters
    const int Z = 10;

    // Plot range
    const double hnuMin = 5 / Constant::ERG_EV;
    const double hnuMax = 15 / Constant::ERG_EV;
    const size_t N = 500;

    ofstream out = IOTools::ofstreamFile("photoelectric/yieldTest.dat");
    for (int i = 0; i < _sizev.size(); i++)
    {
        double a = _sizev[i];
        out << "# a = " << a << '\n';

        // Quantities independent of nu
        double ip_v = ionizationPotential(i, Z);

        double Emin{WD01::eMin(a, Z)};

        // WD01 eq 6
        double hnu_pet = Z >= -1 ? ip_v : ip_v + Emin;

        double hnu = hnuMin;
        const double step = (hnuMax - hnuMin) / N;
        for (size_t n = 0; n < N; n++)
        {
            double hnuDiff = hnu - hnu_pet;
            if (hnuDiff > 0) out << hnu * Constant::ERG_EV << '\t' << photoelectricYield(i, Z, hnuDiff, Emin) << '\n';
            hnu += step;
        }
        out << '\n';
    }
    out.close();
    return 0.0;
}

namespace
{
    // Wavelength grid to use for the tests
    void testSpectrum(double G0, Array& frequencyv, Array& specificIntensityv)
    {
        const double minWav{0.0912 * Constant::UM_CM};  // cutoff at 13.6 eV
        const double maxWav{1000 * Constant::UM_CM};
        const double Tc{3.e4};
        frequencyv = Testing::generateGeometricGridv(200, Constant::LIGHT / maxWav, Constant::LIGHT / minWav);
        specificIntensityv = RadiationFieldTools::generateSpecificIntensityv(frequencyv, Tc, G0);
    }
}  // namespace

void GrainPhotoelectricCalculator::heatingRateTest(double G0, double gasT, double ne) const
{
    Array frequencyv, specificIntensityv;
    testSpectrum(G0, frequencyv, specificIntensityv);

    // Gather environment parameters
    const Spectrum specificIntensity(frequencyv, specificIntensityv);
    Locals env(&specificIntensity, gasT, ne, {-1, 1}, {ne, ne}, {Constant::ELECTRONMASS, Constant::PROTONMASS});

    // File that writes out the absorption efficiency, averaged using the input radiation field as
    // weights.
    ofstream avgQabsOf = IOTools::ofstreamFile("photoelectric/avgQabsInterp.txt");

    // Output file will contain one line for every grain size
    stringstream efficiencyFnSs;
    efficiencyFnSs << "photoelectric/efficiencyG" << setprecision(4) << scientific << G0 << ".dat";
    ofstream efficiencyOf = IOTools::ofstreamFile(efficiencyFnSs.str());

    bool car = true;
    const std::vector<Array> qAbsvv = Testing::qAbsvvForTesting(car, _sizev, frequencyv);

    // For every grain size
    for (size_t i = 0; i < _sizev.size(); i++)
    {
        double a = _sizev[i];
        const Array& Qabsv = qAbsvv[i];

        // Integrate over the radiation field
        Array intensityTimesQabsv = Qabsv * specificIntensityv;
        double intensityQabsIntegral = TemplatedUtils::integrate<double>(frequencyv, intensityTimesQabsv);

        // Calculate and write out the charge distribution and heating efficiency
        ChargeDistribution cd;
        calculateChargeDistribution(i, env, Qabsv, cd);
        stringstream filename;
        filename << "photoelectric/multi-fz/fz_a" << setfill('0') << setw(8) << setprecision(2) << fixed
                 << a / Constant::ANG_CM << ".txt";
        stringstream header;
        header << "# a = " << a << '\n';
        header << "# ne = " << env._ne << '\n';
        header << "# Tgas = " << env._T << '\n';
        cd.plot(filename.str(), header.str());

        double heating = GrainPhotoelectricCalculator::heatingRateA(i, env, Qabsv, cd);

        double totalAbsorbed = Constant::PI * a * a * Constant::FPI * intensityQabsIntegral;
        double efficiency = heating / totalAbsorbed;
        if (!isfinite(efficiency)) cout << "Heating " << heating << " totalabsorbed " << totalAbsorbed << endl;

        efficiencyOf << a / Constant::ANG_CM << '\t' << efficiency << '\n';

        // Calculate and write out the ISRF-averaged absorption efficiency
        double intensityIntegral = TemplatedUtils::integrate<double>(frequencyv, specificIntensityv);
        double avgQabs = intensityQabsIntegral / intensityIntegral;
        avgQabsOf << a / Constant::ANG_CM << '\t' << avgQabs << endl;
    }
    efficiencyOf.close();
    cout << "Wrote " << efficiencyFnSs.str() << endl;
    avgQabsOf.close();
    cout << "Wrote avgQabsInterp.txt" << endl;
    cout << "Charging parameter = " << G0 * sqrt(gasT) / ne << endl;
}

void GrainPhotoelectricCalculator::chargeBalanceTest(double G0, double gasT, double ne, double np) const
{
    Array frequencyv, specificIntensityv;
    testSpectrum(G0, frequencyv, specificIntensityv);
    Spectrum specificIntensity(frequencyv, specificIntensityv);
    Locals env(&specificIntensity, gasT, ne, {-1, 1}, {ne, np}, {Constant::ELECTRONMASS, Constant::PROTONMASS});

    double a = _sizev[0];
    // Qabs for each frequency
    bool car = true;
    Array Qabsv = Testing::qAbsvvForTesting(car, {a}, frequencyv)[0];

    ChargeDistribution cd;
    calculateChargeDistribution(0, env, Qabsv, cd);

    cout << "Zmax = " << cd.zmax() << " Zmin = " << cd.zmin() << '\n';

    stringstream header;
    header << "# a = " << a << endl;
    header << "# G0 = " << G0 << endl;
    header << "# ne = " << ne << endl;
    header << "# Tgas = " << gasT << endl;
    cd.plot("photoelectric/fZ.txt", header.str());
}

double GrainPhotoelectricCalculator::ionizationPotential(int i, int z) const
{
    return WD01::ionizationPotential(_sizev[i], z, _carOrSil);
}

double GrainPhotoelectricCalculator::photoelectricYield(int i, int z, double hnuDiff, double Emin) const
{
    return WD01::yield_cached(_sizev[i], z, hnuDiff, Emin, _carOrSil, _y1Cache[i]);
}

double GrainPhotoelectricCalculator::autoIonizationThreshold(int i) const
{
    return WD01::autoIonizationThreshold(_sizev[i], _carOrSil);
}

double GrainPhotoelectricCalculator::stickingCoefficient(int i, int z, int z_i) const
{
    return WD01::stickingCoefficient_cached(_sizev[i], z, z_i, _carOrSil, _eStickPositiveCache[i],
                                            _eStickNegativeCache[i]);
}
