#ifndef CORE_DEBUGMACROS_H_
#define CORE_DEBUGMACROS_H_

#include <iostream>
#include <ios>

/* #define SANITY */

#ifdef SILENT
#define DEBUG(x)                                                                               \
	do                                                                                     \
	{                                                                                      \
	} while (0)
#else

/* #define PRINT_LEVEL_MATRICES */
/* #define PRINT_CHEMISTRY_MATRICES */
/* #define REPORT_LINE_QUALITY */
/* #define DEBUG_CONTINUUM_DATA */
/* #define ECHO_READIN */

#define DEBUG(x)                                                                               \
	do                                                                                     \
	{                                                                                      \
		std::cout << std::scientific << x;                                            \
	} while (0)
#endif

#endif /* CORE_DEBUGMACROS_H_ */
