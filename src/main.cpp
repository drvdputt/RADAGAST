#include <cstdlib>
#include <string>

#include "IonizationBalance.h"
#include "PhotoelectricHeating.h"
#include "Testing.h"

int main()

{
	//  printf("%f", PhotoelectricHeatingRecipe::yieldFunctionTest());
	//  printf("%f", PhotoelectricHeatingRecipe::chargeBalanceTest());
	//  printf("%f\n", Ionization::testIonization());

	//Testing::testGasSpecies();
	//Testing::testIonizationCrossSection();
	Testing::testPhotoelectricHeating();
	return 0;
}
