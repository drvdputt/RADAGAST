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
	for (double T = 10; T < 100000; T *= 1.1)
	{

		GasInterfaceImpl::Solution s = pimpl->calculateDensities(n, T, specificIntensity, gri, nullptr, kGrainH2);
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
		std::cout << "T = " << T << std::endl;
	}
	outfile.close();
	std::cout << "wrote to " << outfname << std::endl;
}

			
	
	