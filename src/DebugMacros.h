#ifndef _GLOBAL_H_
#define _GLOBAL_H_

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

#endif /* _GLOBAL_H_ */
