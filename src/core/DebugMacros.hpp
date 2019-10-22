#ifndef CORE_DEBUGMACROS_HPP
#define CORE_DEBUGMACROS_HPP

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
#endif // CORE_DEBUGMACROS_HPP
