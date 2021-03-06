#include "SimpleColumnFile.hpp"
#include "IOTools.hpp"
#include <ios>
#include <sstream>

namespace RADAGAST
{
    void InColumnFile::read(int numCols, int reserveLines)
    {
        _columnv.resize(numCols);
        for (auto& v : _columnv) v.reserve(reserveLines);

        std::ifstream ifs = IOTools::ifstreamRepoFile(_fname);
        std::string line;
        while (getline(ifs, line))
        {
            if (line.front() == '#') continue;

            auto iss = std::istringstream(line);
            for (int i = 0; i < numCols; i++)
            {
                double temp;
                iss >> temp;
                if (iss.fail()) Error::runtime("Problem reading column file");

                _columnv[i].emplace_back(temp);
            }
        }

        // In case the given reservation was too large
        for (auto& v : _columnv) v.shrink_to_fit();
    }

    OutColumnFile::OutColumnFile(const std::string& filePath, const std::vector<std::string>& colNamev, int precision)
        : _outFile{IOTools::ofstreamFile(filePath)}, _numCols{colNamev.size()}
    {
        if (precision != -1) _outFile.precision(precision);

        _outFile << std::scientific << '#';
        for (size_t i = 0; i < _numCols; i++) _outFile << i << ' ' << colNamev[i] << ';';
        _outFile << '\n';
    }

    OutColumnFile::~OutColumnFile() { _outFile.close(); }
}
