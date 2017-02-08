#include <cstdlib>
#include <string>

#include "IonizationBalance.h"
#include "PhotoelectricHeating.h"

int main(int argc, char* argv[])

{
    //  printf("%f", Photoelectric::yieldFunctionTest());
    	printf("%f\n", Photoelectric::heatingRateTest("/Users/drvdputt/Testing/ratetest-chargeparams/rateTest_1e4.txt"));
    //  printf("%f", Photoelectric::chargeBalanceTest());
    //  printf("%f\n", Ionization::testIonization());
    return 0;
}
