#ifndef _IOTOOLS_H_
#define _IOTOOLS_H_

#include <fstream>
#include <iostream>
#include <sstream>

/* Just a few shorthands for often repeated statements */

namespace IOTools {
	std::ifstream ifstreamFile(const std::string& file);

	std::istringstream istringstreamNextLine(std::ifstream& ifs);
}



#endif /* _IOTOOLS_H_ */
