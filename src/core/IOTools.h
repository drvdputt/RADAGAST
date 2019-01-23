#ifndef _IOTOOLS_H_
#define _IOTOOLS_H_

#include <fstream>
#include <sstream>
#include <vector>

/* Just a few shorthands for often repeated statements */

namespace IOTools
{
/** Opens a file that resided somewhere in the repo. The macro REPOROOT must be properly set to the
    path of the git directory (without trailing slash, e.g. /home/user/GasModule/git). The correct
    way to provide this option to the compiler is as follows: -DREPOROOT=\""/path/to/git"\". The
    backslashes and quotes need to be there to make sure the preprocessor inserts a string
    literal. */
std::ifstream ifstreamRepoFile(const std::string& pathRelativeToGitDir);

/** Creates an ifstream and opens it using the given filename. An error is thrown if the file could
    not be opened. */
std::ifstream ifstreamFile(const std::string& file);

/** Creates an ofstream and opens it using the given filename. An error is thrown if the file could
    not be opened. */
std::ofstream ofstreamFile(const std::string& file);

/** Reads in a line from a certain ifstream, and immediately creates an istringstream object so that
    it can be parsed. This function moves the ifstream forward one line. */
std::istringstream istringstreamNextLine(std::ifstream& ifs);

std::vector<double> allNumbersFromNextLine(const std::string& line);
} // namespace IOTools

class ColumnFile
{
public:
	ColumnFile(const std::string& filePath, const std::vector<std::string>& colNamev);
	~ColumnFile();
	void writeLine(const std::vector<double>& colValuev);

private:
	std::ofstream _outFile;
	size_t _numCols;
};

#endif /* _IOTOOLS_H_ */
