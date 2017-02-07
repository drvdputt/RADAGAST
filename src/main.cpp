#include "photoelectricHeating.h"
#include "ionizationBalance.h"

#include <cstdlib>

int main(int argc, char* argv[])

{
    //  printf("%f", Photoelectric::yieldFunctionTest());
    printf("%f\n", Photoelectric::heatingRateTest());
    //  printf("%f", Photoelectric::chargeBalanceTest());
    //  printf("%f\n", Ionization::testIonization());
    std::exit(0);
}
