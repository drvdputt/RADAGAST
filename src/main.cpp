#include "Testing.h"

#include <execinfo.h>

#include <unistd.h>
#include <csignal>
#include <iostream>

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
//		Testing::testPhotoelectricHeating();
//		Testing::testIonizationStuff();
//		Testing::testPS64Collisions();
//		Testing::testACollapse();
//		Testing::compareFromFilesvsHardCoded();
		Testing::runFromFilesvsHardCoded();
		Testing::runFullModel();
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
