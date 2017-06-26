#ifndef _FLAGS_H_
#define _FLAGS_H_

#include <string>

#ifndef REPOROOT
#error "Please specify the main GasModule git directory on the compiler command line using -D\"/full/path/to/GasModule/git\""
#endif
const std::string repoRoot = REPOROOT;

#define SANITY

#ifdef SILENT
#define DEBUG(x)                                                                                   \
	do                                                                                         \
	{                                                                                          \
	} while (0)
#else

//#define PRINT_MATRICES
//
//#define REPORT_LINE_QUALITY
//
//#define DEBUG_CONTINUUM_DATA

#define DEBUG(x)                                                                                   \
	do                                                                                         \
	{                                                                                          \
		std::cout << x;                                                                    \
	} while (0)
#endif

#endif /* _FLAGS_H_ */
