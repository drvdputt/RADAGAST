#ifndef _FLAGS_H_
#define _FLAGS_H_

#define SILENT
#ifdef SILENT
	#define DEBUG(x) do {} while (0)
#else
	#define VERBOSE

//	#define PRINT_MATRICES

	#define REPORT_LINE_QUALITY

	#define REPORT_OVERRIDE

	//#define PRINT_CONTINUUM_DATA

	#define DEBUG(x) do {std::cout << x;} while (0)
#endif

#endif /* _FLAGS_H_ */
