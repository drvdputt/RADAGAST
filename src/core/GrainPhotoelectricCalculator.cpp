#include "GrainPhotoelectricCalculator.hpp"
#include "Constants.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"
#include "Functions.hpp"
#include "IOTools.hpp"
#include "Options.hpp"
#include "RadiationFieldTools.hpp"
#include "SpeciesIndex.hpp"
#include "TemplatedUtils.hpp"
#include "Testing.hpp"
#include "WeingartnerDraine2001.hpp"
#include <cmath>

using namespace std;

namespace GasModule
{
    GrainPhotoelectricCalculator::GrainPhotoelectricCalculator(const Array* sizev, const std::vector<Array>* qAbsvv,
                                                               double workFunction, bool carOrSil,
                                                               const Spectrum* meanIntensity)
        : _meanIntensity{meanIntensity},
          _workFunction{workFunction}, _carOrSil{carOrSil}, _sizev{sizev}, _qAbsvv{qAbsvv},
          _integrationWorkspace(meanIntensity->frequencyv().size())
    {
        _y1Cache.resize(_sizev->size());
        _eStickPositiveCache.resize(_sizev->size());
        _eStickNegativeCache.resize(_sizev->size());
        for (size_t i = 0; i < _sizev->size(); i++)
        {
            _y1Cache[i] = WD01::y1((*_sizev)[i]);
            _eStickPositiveCache[i] = WD01::estick_positive((*_sizev)[i]);
            _eStickNegativeCache[i] = WD01::estick_negative((*_sizev)[i]);
        }
    }

    // these are the same for all instances of Locals
    const int GrainPhotoelectricCalculator::Locals::_ie{0};
    const int GrainPhotoelectricCalculator::Locals::_ip{1};
    const int GrainPhotoelectricCalculator::Locals::_iH{2};
    const int GrainPhotoelectricCalculator::Locals::_iH2{3};
    const std::array<int, 4> GrainPhotoelectricCalculator::Locals::_chargev{-1, 1, 0, 0};
    const std::array<double, 4> GrainPhotoelectricCalculator::Locals::_massv{
        Constant::ELECTRONMASS, Constant::PROTONMASS, Constant::HMASS, 2 * Constant::HMASS};

    GrainPhotoelectricCalculator::Locals::Locals(double T, const SpeciesVector& sv) : _T(T)
    {
        // these are different every time a new Locals is created
        _densityv = {sv.ne(), sv.np(), sv.nH(), sv.nH2()};
    }

    int GrainPhotoelectricCalculator::minimumCharge(int i) const
    {
        double Uait = autoIonizationThreshold(i);
        return WD01::minimumCharge((*_sizev)[i], Uait);
    }

    void GrainPhotoelectricCalculator::calculateChargeDistribution(int i, Locals& env, ChargeDistribution& cd)
    {
        double a = (*_sizev)[i];

        // Express a in angstroms
        double aA = a / Constant::ANGSTROM;

        // Highest possible energy of a photon
        const Array& frequencyv = _meanIntensity->frequencyv();
        double hnumax = Constant::PLANCK * frequencyv[frequencyv.size() - 1];

        // The maximum charge is one more than the highest charge which still allows ionization by
        // photons of hnumax.
        int resultZmax = floor(((hnumax - _workFunction) / (14.4 * Constant::EV) * aA + .5 - .3 / aA) / (1 + .3 / aA));

        // The minimum charge is the most negative charge for which autoionization does not occur
        int resultZmin = minimumCharge(i);

        // The few cases I've seen this happen, Zmin is always 0 and Zmax is always -1. Since Zmax
        // only depends on the maximum photon energy, it should be fine to increase it up to Zmin
        // in these edge cases.
        if (resultZmax < resultZmin) resultZmax = resultZmin;

        // These two functions determine the up and down rates for the detailed balance

        // The rate at which the grain moves out of charge Z, in the positive direction.
        // Contributions by photoelectric ejection of electron, and collisions with positive
        // particles. Collisions with multiply charged particles are treated as if they can only
        // change the charge in steps of 1. This is technically incorrect, but the charging rate
        // due to positive ions is very small anyway.
        auto chargeUpRate = [&](int z) {
            double Jtotal{0};
            Jtotal += emissionRate(i, z);
            for (size_t j = 0; j < env._chargev.size(); j++)
            {
                if (env._chargev[j] > 0)
                    Jtotal += collisionalChargingRate(i, env._T, z, env._chargev[j], env._massv[j], env._densityv[j]);
            }
            return Jtotal;
        };
        // The rate at which the grain moves out of charge Z, in the negative direction. Sums
        // collisions over negative particles.
        auto chargeDownRate = [&](int z) {
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
    }

    void GrainPhotoelectricCalculator::getPET_PDT_Emin(int i, int Z, double& pet, double& pdt, double& Emin) const
    {
        Emin = WD01::eMin((*_sizev)[i], Z);
        double ip_v = ionizationPotential(i, Z);
        // WD01 eq 6
        pet = Z >= -1 ? ip_v : ip_v + Emin;
        // Photodetachment, WD01 eq 18 (using EA(Z + 1, a) = IP(Z, a) if Z < 0)
        pdt = ip_v + Emin;
    }

    double GrainPhotoelectricCalculator::photoelectricIntegrationLoop(int i, double nuPET, const Array& fNu)
    {
        const Array& frequencyv = _meanIntensity->frequencyv();
        const Array& meanIntensityv = _meanIntensity->valuev();
        const Array& Qabsv = _qAbsvv->at(i);

        // yield is zero by definition below photoelectric threshold
        int j = frequencyv.size() - 1;
        while (j >= 0 && frequencyv[j] > nuPET)
        {
            double hnu = Constant::PLANCK * frequencyv[j];
            _integrationWorkspace[j] = Qabsv[j] * meanIntensityv[j] / hnu * fNu[j];
            j--;
        }

        if (j >= frequencyv.size() - 1)
            return 0;
        else
        {
            // don't forget to increment i!

            // <function unit> sr-1 cm-2 s-1
            double integral =
                TemplatedUtils::integrate<double>(frequencyv, _integrationWorkspace, j + 1, frequencyv.size() - 1);
            // <function unit> cm-2 s-1
            return Constant::FPI * integral;
        }
    }

    double GrainPhotoelectricCalculator::photodetachmentIntegrationLoop(int Z, double pdt,
                                                                        const double* calcEnergyWithThisEmin)
    {
        const Array& frequencyv = _meanIntensity->frequencyv();
        const Array& meanIntensityv = _meanIntensity->valuev();

        // no effect below photodetachment threshold
        double nu_pdt = pdt / Constant::PLANCK;
        int j = frequencyv.size() - 1;
        while (j >= 0 && frequencyv[j] > nu_pdt)
        {
            double hnu = Constant::PLANCK * frequencyv[j];
            double hnuDiff = hnu - pdt;
            // <function unit> / time / angle
            _integrationWorkspace[j] = WD01::sigmaPDT(Z, hnuDiff) * meanIntensityv[j] / hnu;
            if (calcEnergyWithThisEmin) _integrationWorkspace[j] *= hnuDiff + *calcEnergyWithThisEmin;
            j--;
        }

        if (j >= frequencyv.size() - 1)
            return 0;
        else
        {
            // <erg optional> s-1 sr-1
            double integral =
                TemplatedUtils::integrate<double>(frequencyv, _integrationWorkspace, j + 1, frequencyv.size() - 1);
            // <erg optional> s-1
            return Constant::FPI * integral;
        }
    }

    double GrainPhotoelectricCalculator::heatingRateAZ(int i, int Z)
    {
        auto cacheEntry = _heatingRateCache.find({i, Z});
        if (cacheEntry != _heatingRateCache.end()) return cacheEntry->second;

        double a = (*_sizev)[i];
        const double e2_a = Constant::ESQUARE / a;

        // It's cheaper to calculate these together (and safer, less code duplication)
        double pet, pdt, Emin;
        getPET_PDT_Emin(i, Z, pet, pdt, Emin);
        double nuPET = pet / Constant::PLANCK;

        // WD01 text between eq 10 and 11
        double Elow = Z < 0 ? Emin : -(Z + 1) * e2_a;

        // use cached yield if available
        auto yieldCacheEntry = _yieldCache.find({i, Z});

        const Array& frequencyv = _meanIntensity->frequencyv();
        Array yieldTimesAverageEnergyv(frequencyv.size());
        for (int j = frequencyv.size() - 1; j >= 0 && frequencyv[j] > nuPET; j--)
        {
            double hnuDiff = Constant::PLANCK * frequencyv[j] - pet;
            double Y = yieldCacheEntry != _yieldCache.cend() ? yieldCacheEntry->second[j]
                                                             : photoelectricYield(i, Z, hnuDiff, Emin);
            double Ehigh = Z < 0 ? hnuDiff + Emin : hnuDiff;
            // The integral over the electron energy distribution (integral E f(E) dE), over the
            // energy range for which electrons can escape
            double IntE = WD01::energyIntegral(Elow, Ehigh, Emin);
            // Divide by (integral f(E) dE) over the same range, which normalizes the above value
            // --> IntE / y2 gives an average energy
            double y2 = WD01::escapingFraction(Z, Elow, Ehigh);
            yieldTimesAverageEnergyv[j] = Y * IntE / y2;
        }
        double heatingRatePE = Constant::PI * a * a * photoelectricIntegrationLoop(i, nuPET, yieldTimesAverageEnergyv);

        double heatingRatePD = 0;
        if (Z < 0) heatingRatePD = photodetachmentIntegrationLoop(Z, pdt, &Emin);

        double result = heatingRatePE + heatingRatePD;
        _heatingRateCache[{i, Z}] = result;
        return result;
    }

    double GrainPhotoelectricCalculator::heatingRateA(int i, const ChargeDistribution& cd)
    {
        return cd.sumOverCharge([&](int z) { return heatingRateAZ(i, z); });
    }

    double GrainPhotoelectricCalculator::emissionRate(int i, int Z)
    {
        auto cacheEntry = _emissionRateCache.find({i, Z});
        if (cacheEntry != _emissionRateCache.end()) return cacheEntry->second;

        double a = (*_sizev)[i];
        double pet, pdt, Emin;
        getPET_PDT_Emin(i, Z, pet, pdt, Emin);
        double nuPET = pet / Constant::PLANCK;

        const Array& frequencyv = _meanIntensity->frequencyv();
        // First check if we have yieldv in cache for this size and charge. If not, calculate it
        // (remember that [] inserts an empty Array if it did not exist yet
        Array& yieldv = _yieldCache[{i, Z}];
        if (yieldv.size() != frequencyv.size()) yieldv.resize(frequencyv.size());
        for (int j = frequencyv.size() - 1; j >= 0 && frequencyv[j] > nuPET; j--)
        {
            double hnuDiff = Constant::PLANCK * frequencyv[j] - pet;
            yieldv[j] = photoelectricYield(i, Z, hnuDiff, Emin);
        }
        double emissionRatePE = Constant::PI * a * a * photoelectricIntegrationLoop(i, nuPET, yieldv);
        double emissionRatePD = 0;
        if (Z < 0) emissionRatePD = photodetachmentIntegrationLoop(Z, pdt, nullptr);

        double result = emissionRatePE + emissionRatePD;
        _emissionRateCache[{i, Z}] = result;
        return result;
    }

    double GrainPhotoelectricCalculator::collisionalChargingRate(int i, double gasT, int Z, int particleCharge,
                                                                 double particleMass, double particleDensity) const
    {
        double a = (*_sizev)[i];
        double kT = Constant::BOLTZMAN * gasT;

        // WD eq 26: akT / q^2 = akT / e^2 / z^2
        double tau = a * kT / Constant::ESQUARE / particleCharge / particleCharge;
        // Ze / q = Z / z
        double nu = static_cast<double>(Z) / static_cast<double>(particleCharge);

        double Jtilde;
        if (nu < 0)
        {
            Jtilde = (1. - nu / tau) * (1. + sqrt(2. / (tau - 2. * nu)));
        }
        else if (nu > 0)
        {
            double toSquare = 1. + 1. / sqrt(4. * tau + 3. * nu);
            Jtilde = toSquare * toSquare * exp(-WD01::thetaNu(nu) / tau);
        }
        else
        {
            Jtilde = 1. + sqrt(Constant::PI / 2. / tau);
        }
        return particleDensity * stickingCoefficient(i, Z, particleCharge) * sqrt(8. * kT * Constant::PI / particleMass)
               * a * a * Jtilde;
    }

    double GrainPhotoelectricCalculator::recombinationCoolingRate(int i, const Locals& env,
                                                                  const ChargeDistribution& cd) const
    {
        double a = (*_sizev)[i];

        // Calculates WD01 equation 42
        double kT = Constant::BOLTZMAN * env._T;
        double eightkT3DivPi = 8 * kT * kT * kT / Constant::PI;

        // For every charged collision partner (effectively electrons and protons), add the
        // contributions for each possible grain charge.
        double particleSum = 0;
        for (size_t j = 0; j < env._chargev.size(); j++)
        {
            int zParticle = env._chargev[j];
            if (zParticle)
            {
                // tau = akT / q^2 (WD01 eq 26)
                double tau = a * kT / (zParticle * zParticle * Constant::ESQUARE);
                double Zsum = cd.sumOverCharge([&](int zGrain) {
                    double ksi = zGrain / static_cast<double>(zParticle);
                    return stickingCoefficient(i, zGrain, zParticle) * WD01::lambdaTilde(tau, ksi);
                });
                particleSum += env._densityv[j] * sqrt(eightkT3DivPi / env._massv[j]) * Zsum;
            }
        }

        // The second term of equation 42: autoionization of grains with the most negative charge
        // inhibits the cooling of the gas. EA(Zmin) = IP(Zmin-1) because IP(Z) = EA(Z+1)
        double secondTerm = 0;
        // This term is only included when the population of the maximally negative grain charge
        // minimumCharge is significant. If it is not siginicant, then fZ will not cover
        // minimumCharge, (and Zmin > minimumCharge).
        int zmin = cd.zmin();
        if (zmin == minimumCharge(i))
            secondTerm = cd.value(zmin)
                         * collisionalChargingRate(i, env._T, zmin, -1, Constant::ELECTRONMASS, env._densityv[env._ie])
                         * ionizationPotential(i, zmin - 1);

        return Constant::PI * a * a * particleSum * kT + secondTerm;
    }

    double GrainPhotoelectricCalculator::gasGrainCollisionCooling(int i, const Locals& env,
                                                                  const ChargeDistribution& cd, double Tgrain,
                                                                  bool forGrain) const
    {
        double a = (*_sizev)[i];
        double kT = env._T * Constant::BOLTZMAN;
        double kTgrain = Tgrain * Constant::BOLTZMAN;

        double lambdaG = cd.sumOverCharge([&](int zGrain) {
            double lambdaG_for_this_z = 0;
            // Not entirely sure if this is the right potential
            double Ug = ionizationPotential(i, zGrain);
            double Vg = sqrt(Constant::ESQUARE) * Ug;
            for (size_t j = 0; j < env._chargev.size(); j++)
            {
                // Dimensionless
                int zParticle = env._chargev[j];
                double ZVg = zParticle * Vg;
                double psi = ZVg / kT;
                double eta = psi <= 0 ? 1 - psi : exp(-psi);
                double ksi = psi <= 0 ? 1 - psi / 2 : (1 + psi / 2) * exp(-psi);
                // incoming collision
                double total = 2 * kT * ksi;

                // energy carried away by evaporation of the neutralized/thermalized version of the
                // colliding particle (does not count for electrons, see text below eq 29)
                if (j != env._ie) total -= 2 * kTgrain * eta;

                if (forGrain)
                {
                    // A charged particle will slow down or speed up before it hits a charged grain
                    total -= ZVg * eta;

                    // when protons charge the grain, this means they recombine on the surface
                    // (not included in chemical network yet though). The recombination energy
                    // needs to go somewhere, so the grain gets extra heat. I'm pretty sure
                    // using 1 rydberg (= 13.6 eV = ionization potential) is fine here.
                    if (j == env._ip) total += Constant::RYDBERG * eta;
                }

                double S = 0.;
                if (zParticle)
                {
                    // For charged particles, use the same sticking coefficient that was used
                    // for the charging rates (for consistency).
                    S = stickingCoefficient(i, zGrain, zParticle);
                }
                else
                {
                    // For neutral particles, S needs to be an accomodation coefficient instead.
                    // Baldwin states that it should be mM / (m^2 + M^2) and references Draine
                    // (1978), with M the mass of a typical atom in the grain.
                    double m = env._massv[j];
                    double M = WD01::atomMass(_carOrSil);
                    S = m * M / (m * m + M * M);
                }

                // cm s-1
                double vbar = Functions::meanThermalVelocity(env._T, env._massv[j]);

                // cm-3 * cm s-1 * erg = cm-2 s-1 erg
                lambdaG_for_this_z += env._densityv[j] * vbar * S * total;
            }
            return lambdaG_for_this_z;
        });
        // Remember that 1991-Baldwin gives the rate 'per unit projected area'
        return Constant::PI * a * a * lambdaG;
    }

    double GrainPhotoelectricCalculator::ionizationPotential(int i, int z) const
    {
        return WD01::ionizationPotential((*_sizev)[i], z, _carOrSil);
    }

    double GrainPhotoelectricCalculator::photoelectricYield(int i, int z, double hnuDiff, double Emin) const
    {
        return WD01::yield_cached((*_sizev)[i], z, hnuDiff, Emin, _carOrSil, _y1Cache[i]);
    }

    double GrainPhotoelectricCalculator::autoIonizationThreshold(int i) const
    {
        return WD01::autoIonizationThreshold((*_sizev)[i], _carOrSil);
    }

    double GrainPhotoelectricCalculator::stickingCoefficient(int i, int z, int z_i) const
    {
        return WD01::stickingCoefficient_cached((*_sizev)[i], z, z_i, _carOrSil, _eStickPositiveCache[i],
                                                _eStickNegativeCache[i]);
    }
}
