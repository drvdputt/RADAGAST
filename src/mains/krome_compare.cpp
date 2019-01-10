#include "GasInterface.h"
#include "GasInterfaceImpl.h"
#include "SpecialFunctions.h"
#include "SpeciesIndex.h"
#include "Spectrum.h"
#include "Testing.h"

#include <fstream>

int main()
{
	// Radiation field
	int nPhotoBins = 20;
	double low_eV = 1;
	double high_eV = 24;
	double Tc = 4000;
	// e = h nu ; n = e / h
	double lowFreq = low_eV / Constant::ERG_EV / Constant::PLANCK;
	double highFreq = high_eV / Constant::ERG_EV / Constant::PLANCK;
	Array frequencyv = Testing::generateGeometricGridv(nPhotoBins, lowFreq, highFreq);
	Array intensityv(frequencyv.size());
	for (int i = 0; i < frequencyv.size(); i++)
		intensityv[i] = SpecialFunctions::planck(frequencyv[i], Tc);
	Spectrum specificIntensity(frequencyv, intensityv);

	// Write out the radiation field in eV / cm2 / sr / hz
	Array ev = Constant::PLANCK * frequencyv * Constant::ERG_EV;
	Array iv = intensityv * Constant::ERG_EV;
	std::ofstream radf("photobinj.dat");
	radf << "# Emid (eV)\tJ (eV / cm2 / sr)\n";
	for (int i = 0; i < nPhotoBins; i++)
	{
		double j = iv[i];
		radf << ev[i] << '\t' << j << '\n';
	}
	radf.close();

	// Set a fixed H2 formation rate
	double kGrainH2 = 1e-12;

	// H nuclei density
	double n = 1e4;

	// Load a model, and grab the underlying implementation
	GasModule::GasInterface gasInterface(frequencyv, frequencyv, frequencyv, "", "5 0");
	auto pimpl = gasInterface.pimpl();

	// No grains, we will manually override the h2 formation rate instead
	GasModule::GrainInterface gri{};

	std::string outfname = "equilibrium_densities";
	std::ofstream outfile(outfname);
	outfile << "# T e H H2 H+ heat cool";
	GasInterfaceImpl::Solution s;
	GasInterfaceImpl::Solution* sp = nullptr;
	for (double T = 1000; T < 100000; T *= 1.05)
	{
		s = pimpl->calculateDensities(n, T, specificIntensity, gri, sp, kGrainH2);
		// Uncomment this to re-use the previous solution as an initial guess
		// sp = &s;
		double heat = pimpl->heating(s, gri);
		double cool = pimpl->cooling(s);
		outfile << T;
		for (const std::string& name : {"e-", "H", "H2", "H+"})
		{
			int i = SpeciesIndex::index(name);
			outfile << " " << s.speciesNv[i];
		}
		outfile << " " << heat << " " << cool;
		outfile << '\n';
		std::cout << "T = " << T << '\n';
	}
	for (double T = 100000; T > 10; T /= 1.05)
	{
		s = pimpl->calculateDensities(n, T, specificIntensity, gri, sp, kGrainH2);
		sp = &s;
		double heat = pimpl->heating(s, gri);
		double cool = pimpl->cooling(s);
		outfile << T;
		for (const std::string& name : {"e-", "H", "H2", "H+"})
		{
			int i = SpeciesIndex::index(name);
			outfile << " " << s.speciesNv[i];
		}
		outfile << " " << heat << " " << cool;
		outfile << '\n';
		std::cout << "T = " << T << '\n';
	}

	outfile.close();
	std::cout << "wrote to " << outfname << '\n';
}
