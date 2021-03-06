#include "IOTools.hpp"
#include "DebugMacros.hpp"
#include "Error.hpp"

#ifndef REPOROOT
#    error \
        "Please specify the main RADAGAST git directory on the compiler command line using -D\"/full/path/to/RADAGAST/git\""
#endif

using namespace std;

namespace
{
    const string repoRoot = REPOROOT;
}

namespace RADAGAST
{
    ifstream IOTools::ifstreamRepoFile(const string& pathRelativeToGitDir)
    {
        return ifstreamFile(repoRoot + "/" + pathRelativeToGitDir);
    }

    ifstream IOTools::ifstreamFile(const string& file)
    {
        ifstream input(file);
        if (!input) Error::runtime("Input file " + file + "not found.");
        DEBUG("Opened file (read) " << file << '\n');
        return input;
    }

    ofstream IOTools::ofstreamFile(const string& file)
    {
        ofstream output(file);
        if (!output.is_open()) Error::runtime("Output file " + file + " could not be opened.");
        DEBUG("Opened file (write) " << file << '\n');
        return output;
    }

    istringstream IOTools::istringstreamNextLine(ifstream& ifs)
    {
        string line;
        getline(ifs, line);
        return istringstream(line);
    }

    vector<double> IOTools::allNumbersFromNextLine(const string& line)
    {
        auto ss = istringstream(line);
        double i;
        vector<double> nv;
        while (ss >> i) nv.emplace_back(i);
        return nv;
    }
}
