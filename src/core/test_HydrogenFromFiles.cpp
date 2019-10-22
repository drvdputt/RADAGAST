#include "DoctestUtils.hpp"

#include "EigenAliases.hpp"
#include "GasStruct.hpp"
#include "HFromFiles.hpp"
#include "HHardCoded.hpp"
#include "SpeciesIndex.hpp"

TEST_CASE("Test correctness of collapsed A-matrix for hydrogen")
{
	/** Writes out the A-coefficients of a fully collapsed H-model, so that they can be compared to the
	    NIST values between different n. */
	HFromFiles hff(0);
	EMatrix avv = hff.avv().array();
	// cout << hff.avv() << endl;
	// cout << "Compare this with values directly from NIST below:" << endl;
	EMatrix nistA(5, 5);
	// clang-format off
	nistA << 0, 0, 0, 0, 0,
		4.6986e+08, 0, 0, 0, 0,
		5.5751e+07, 4.4101e+07, 0, 0, 0,
		1.2785e+07, 8.4193e+06, 8.9860e+06, 0, 0,
		4.1250e+06, 2.5304e+06, 2.2008e+06, 2.6993e+06, 0;
	// clang-format on
	// cout << nistA << endl;
	// cout << "The element-wise relative difference is " << endl;
	EMatrix relDiff = (avv - nistA).array() / nistA.array();
	DoctestUtils::compareMatrices(avv, nistA, 0.002);
}

TEST_CASE("Perform a direct comparison of the two LevelDataProvider classes HydrogenHardcoded "
          "and HydrogenFromFiles.")
{
	HydrogenHardcoded hhc;
	HFromFiles hff(2);

	// auto hc_vs_ff = [&](auto hc_thing, auto ff_thing) {
	// 	cout << "Hardcoded" << endl;
	// 	cout << hc_thing << endl;
	// 	cout << "From files" << endl;
	// 	cout << ff_thing << endl;
	// };

	assert(hhc.numLv() == hff.numLv());

	// cout << "Energy levels:" << endl;
	EVector evhc = hhc.ev();
	EVector evff = hff.ev();
	// hc_vs_ff(evhc, evff);
	DoctestUtils::compareMatrices(evhc, evff, 0.01);

	assert(hhc.gv() == hff.gv());

	// cout << "A coefficients:" << endl;
	EMatrix avvhc = hhc.avv();
	EMatrix avvff = hff.avv();
	// hc_vs_ff(avvhc, avvff);
	DoctestUtils::compareMatrices(avvhc, avvff, 0.01);

	// cout << "Extra A:" << endl;
	EMatrix eavvhc = hhc.extraAvv();
	EMatrix eavvff = hff.extraAvv();
	// hc_vs_ff(eavvhc, eavvff);
	DoctestUtils::compareMatrices(eavvhc, eavvff, 0.01);

	double T = 1e4;
	double ne = 1e4;
	double np = 1e4;
	EVector speciesNv{EVector::Zero(SpeciesIndex::size())};
	speciesNv(SpeciesIndex::inp()) = np;
	speciesNv(SpeciesIndex::ine()) = ne;
	GasStruct gas(T, speciesNv);

	// cout << "Collisions:" << endl;
	EMatrix cvvhc = hhc.cvv(gas);
	EMatrix cvvff = hff.cvv(gas);
	// hc_vs_ff(cvvhc, cvvff);
}
