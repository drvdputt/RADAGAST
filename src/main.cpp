#include <cstdlib>
#include <string>

#include "IonizationBalance.h"
#include "PhotoelectricHeating.h"
#include "Testing.h"

int main()

{
	//  printf("%f", Photoelectric::yieldFunctionTest());
	//Photoelectric::heatingRateTest("/Users/drvdputt/GasModule/run/rateTest.dat");
	//  printf("%f", Photoelectric::chargeBalanceTest());
	//  printf("%f\n", Ionization::testIonization());

	Testing::testGasSpecies();

	return 0;
}
