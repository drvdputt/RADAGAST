#include "Functions.hpp"
#include "GasInterface.hpp"
#include "GasInterfaceImpl.hpp"
#include "SpeciesIndex.hpp"
#include "Spectrum.hpp"
#include "Testing.hpp"
#include <fstream>

using namespace RADAGAST;

int main()
{
    SpeciesIndex spindex(SpeciesIndex::e_p_H_H2);

    // Radiation field
    int nPhotoBins = 20;
    double low_eV = 1;
    double high_eV = 24;
    double Tc = 2950;
    // e = h nu ; n = e / h
    double lowFreq = low_eV * Constant::EV / Constant::PLANCK;
    double highFreq = high_eV * Constant::EV / Constant::PLANCK;
    Array frequencyv = Testing::generateGeometricGridv(nPhotoBins, lowFreq, highFreq);
    Array intensityv(frequencyv.size());
    for (int i = 0; i < frequencyv.size(); i++) intensityv[i] = Functions::planck(frequencyv[i], Tc);
    Spectrum meanIntensity(frequencyv, intensityv);

    // Write out the radiation field in eV / cm2 / sr / hz
    Array ev = Constant::PLANCK * frequencyv / Constant::EV;
    Array iv = intensityv / Constant::EV;
    std::ofstream radf("photobinj.dat");
    radf << "# Emid (eV)\tJ (eV / cm2 / sr)\n";
    for (int i = 0; i < nPhotoBins; i++)
    {
        double j = iv[i];
        radf << ev[i] << '\t' << j << '\n';
    }
    radf.close();

    // Set a fixed H2 formation rate
    double kGrainH2 = 2e-15;

    // H nuclei density
    double n = 1e4;

    // Load a model, and grab the underlying implementation
    RADAGAST::GasInterface gasInterface(frequencyv, frequencyv, frequencyv);

    // No grains, we will manually override the h2 formation rate instead
    RADAGAST::GrainInterface gri{};

    std::string outfname = "equilibrium_densities";
    std::ofstream outfile(outfname);
    outfile << "# T e H H2 H+ heat cool";
    for (double T = 1000; T < 100000; T *= 1.05)
    {
        GasSolution s = gasInterface.solveDensities(n, T, meanIntensity, 1., gri, kGrainH2);
        // Uncomment this to re-use the previous solution as an initial guess
        // sp = &s;
        double heat = s.heating();
        double cool = s.cooling();
        outfile << T;
        for (const std::string& name : {"e-", "H", "H2", "H+"})
        {
            outfile << " " << s.speciesVector().nSpecies(name);
        }
        outfile << " " << heat << " " << cool;
        outfile << '\n';
        std::cout << "T = " << T << '\n';
    }
    for (double T = 100000; T > 10; T /= 1.05)
    {
        GasSolution s = gasInterface.solveDensities(n, T, meanIntensity, 1., gri, kGrainH2);
        // sp = &s;
        double heat = s.heating();
        double cool = s.cooling();
        outfile << T;
        for (const std::string& name : {"e-", "H", "H2", "H+"})
        {
            outfile << " " << s.speciesVector().nSpecies(name);
        }
        outfile << " " << heat << " " << cool;
        outfile << '\n';
        std::cout << "T = " << T << '\n';
    }

    outfile.close();
    std::cout << "wrote to " << outfname << '\n';
}
