#include "Testing.h"

#include <csignal>
#include <execinfo.h>
#include <iostream>
#include <unistd.h>

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
	signal(SIGSEGV, handler);
	try
	{
		// These output graphable data which can be used to check correctness.
		// Testing::plotPhotoelectricHeating();
		// Testing::plotIonizationStuff();
		// Testing::plotPS64Collisions();
		// Testing::plotInterpolationTests();

		// These do a run, with writeout, but the results are not directly comparable to
		// something in the way that the ones above work.
		// Testing::runFromFilesvsHardCoded();
		// Testing::runFullModel();
		Testing::runWithDust(true);
		// Testing::runH2(false);
	}
	catch (const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}
	catch (char const* str)
	{
		std::cerr << str << std::endl;
	}
	return 0;
}
