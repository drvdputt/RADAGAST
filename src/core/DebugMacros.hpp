#ifndef CORE_DEBUGMACROS_H_
#define CORE_DEBUGMACROS_H_

#include <iostream>
#include <ios>

#ifdef SILENT
#define DEBUG(x)                                                                               \
	do                                                                                     \
	{                                                                                      \
	} while (0)
#else
#define DEBUG(x)                                                                               \
	do                                                                                     \
	{                                                                                      \
		std::cout << std::scientific << x;                                            \
	} while (0)
#endif
#endif /* CORE_DEBUGMACROS_H_ */
