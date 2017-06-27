#ifndef _IOTOOLS_H_
#define _IOTOOLS_H_

#include <fstream>
#include <iostream>
#include <sstream>

/* Just a few shorthands for often repeated statements */

namespace IOTools
{
/* Creates an ifstream and opens it using the given filename. An error is thrown if the file could
   not be opened. */
std::ifstream ifstreamFile(const std::string& file);

/* Creates an ofstream and opens it using the given filename. An error is thrown if the file could
   not be opened. */
std::ofstream ofstreamFile(const std::string& file);

/* Reads in a line from a certain ifstream, and immediately creates an istringstream object so that
   it can be parsed. This function moves the ifstream forward one line. */
std::istringstream istringstreamNextLine(std::ifstream& ifs);
}

#endif /* _IOTOOLS_H_ */
