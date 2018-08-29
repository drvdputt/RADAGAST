/** This file contains a main which will run all tests. Or something. */

#include "LineProfile.h"
#include "Testing.h"

int main()
{
	test_addToSpectrum();
	test_integrateSpectrum();

	Testing::testChemistry();
	Testing::testACollapse();
	Testing::testFromFilesvsHardCoded();
}
