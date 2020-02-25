#ifndef CORE_IOTOOLS_HPP
#define CORE_IOTOOLS_HPP

#include <fstream>
#include <sstream>
#include <vector>

/** Some shorthand for opening file handles. */
namespace IOTools
{
    /** Opens a file that resided somewhere in the repo. The macro REPOROOT must be properly set to
        the path of the git directory (without trailing slash, e.g. /home/user/GasModule/git). The
        correct way to provide this option to the compiler is as follows:
        -DREPOROOT=\""/path/to/git"\". The backslashes and quotes need to be there to make sure the
        preprocessor inserts a string literal. */
    std::ifstream ifstreamRepoFile(const std::string& pathRelativeToGitDir);

    /** Creates an ifstream and opens it using the given filename. An error is thrown if the file
        could not be opened. */
    std::ifstream ifstreamFile(const std::string& file);

    /** Creates an ofstream and opens it using the given filename. An error is thrown if the file
        could not be opened. */
    std::ofstream ofstreamFile(const std::string& file);

    /** Reads in a line from a certain ifstream, and immediately creates an istringstream object so
        that it can be parsed. This function moves the ifstream forward one line. */
    std::istringstream istringstreamNextLine(std::ifstream& ifs);

    std::vector<double> allNumbersFromNextLine(const std::string& line);
}  // namespace IOTools

#endif  // CORE_IOTOOLS_HPP
