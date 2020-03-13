#include "Testing.hpp"
#include <execinfo.h>
#include <unistd.h>
#include <csignal>
#include <iostream>
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

int main(int argc, char** argv)
{
    if (!(argc == 1 || argc == 4))
    {
        std::cerr << "wrong number of command line arguments\nEither use no arguments "
                     "or specify nH, Tc and lumSol\n ";
        exit(1);
    }

    signal(SIGSEGV, handler);
    try
    {
        // // These output graphable data which can be used to check correctness.
        // Testing::plotPhotoelectricHeating();
        // Testing::plotIonizationStuff();
        // Testing::plotPS64Collisions();
        // Testing::plotInterpolationTests();
        // Testing::plotHeatingCurve_main();

        // These do a run, with writeout, but the results are not directly comparable to
        // something in the way that the ones above work.
        // Testing::runFromFilesvsHardCoded();
        // Testing::runFullModel();
        // Testing::runH2(true);
        if (argc == 1)
            GasModule::Testing::runMRNDust(true);
        else if (argc == 4)
        {
            double nH = std::stod(argv[1]);
            double Tc = std::stod(argv[2]);
            double lumSol = std::stod(argv[3]);
            GasModule::Testing::runMRNDust(true, nH, Tc, lumSol, false);
        }
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
