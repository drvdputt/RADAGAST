#ifndef _FLAGS_H_
#define _FLAGS_H_

//#define SILENT

//#define PRINT_CONTINUUM_DATA

#define VERBOSE

#define PRINT_MATRICES

#define REPORT_LINE_QUALITY

#ifdef SILENT
	#define DEBUG(x)
#else
	#define DEBUG(x) do {std::cout << x;} while (0)
#endif

#endif /* _FLAGS_H_ */
