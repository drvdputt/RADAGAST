#include "DoctestUtils.hpp"
#include "Ionization.hpp"
#include "RadiationFieldTools.hpp"
#include "SimpleHChemistry.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"

using namespace RADAGAST;

TEST_CASE("SimpleHChemistry: compare exact solution of ionization")
{
    const double T = 10000;
    Array frequencyv = Testing::generateGeometricGridv(200, 1e11, 1e16);
    Array meanIntensityv = RadiationFieldTools::generateSpecificIntensityv(frequencyv, 25000, 10);
    Spectrum meanIntensity(frequencyv, meanIntensityv);

    SimpleHChemistry chemistry{};

    // no h2
    double kform = 0;
    double kdiss = 0;
    EVector kv = chemistry.rateCoeffv(T, meanIntensity, kdiss, kform);
    std::stringstream ss;
    ss << "Rate coeff: ionization, recombination, dissociation\n" << kv << '\n';
    std::string ratesMessage = ss.str();
    CAPTURE(ratesMessage);

    SpeciesVector sv0(&chemistry.speciesIndex());
    sv0.setDensities(chemistry.speciesIndex().linearCombination(SpeciesIndex::e_p_H_H2, {50, 50, 50, 0}));
    EVector nv = chemistry.solveBalance(kv, sv0.speciesNv());
    SpeciesVector sv(&chemistry.speciesIndex());
    sv.setDensities(nv);
    double f_network = sv.np() / (sv.np() + sv.nH());
    double f_exact = Ionization::solveBalance(sv0.nH() + sv0.np(), T, meanIntensity);
    CHECK_MESSAGE(TemplatedUtils::equalWithinTolerance(f_exact, f_network, 1e-6),
                  "exact ionization rate = " << f_exact << " while the one from chemical network is " << f_network
                                             << " (species vector is " << nv << ")");
}
