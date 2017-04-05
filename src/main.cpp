#include "PhotoelectricHeating.h"
#include "Testing.h"

#include <execinfo.h>
#include <unistd.h>

#include <cstdlib>
#include <csignal>
#include <string>

void handler(int sig)
{
	void* array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}

int main()
{
	//signal(SIGSEGV, handler);
	//  printf("%f", PhotoelectricHeatingRecipe::yieldFunctionTest());
	//  printf("%f", PhotoelectricHeatingRecipe::chargeBalanceTest());
	//  printf("%f\n", Ionization::testIonization());

	Testing::testGasInterfaceImpl();
	// Testing::testIonizationCrossSection();
	// Testing::testPhotoelectricHeating();
	return 0;
}
