#ifndef CORE_SIMPLECOLUMNFILE_HPP
#define CORE_SIMPLECOLUMNFILE_HPP

#include "Error.hpp"
#include <fstream>
#include <string>
#include <vector>

namespace RADAGAST
{
    /** Class to read in data file with tab or space separated columns. */
    class InColumnFile
    {
    public:
        /** Create a new instance, setting the file name. No file handle is opened yet */
        InColumnFile(const std::string& fname) : _fname{fname} {}

        /** Store the data from the first numCols columns. For speed, a guess for the number of
            lines can be given. */
        void read(int numCols, int reserveLines = 0);

        /** Get the data from one of the columns */
        const std::vector<double>& column(int i) { return _columnv[i]; }

    private:
        std::string _fname;
        std::vector<std::vector<double>> _columnv;
    };

    class OutColumnFile
    {
    public:
        OutColumnFile(const std::string& filePath, const std::vector<std::string>& colNamev, int precision = -1);
        ~OutColumnFile();

        /** Writes one line to the file, separated by spaces. The argument should be of a typical
            container type, preferably of doubles (supporting iterators and size). */
        template<typename T> void writeLine(const T& colValuev);

    private:
        std::ofstream _outFile;
        size_t _numCols;
    };

    template<typename T> void OutColumnFile::writeLine(const T& colValuev)
    {
        Error::equalCheck("numCols and num values in line", _numCols, colValuev.size());
        auto it = begin(colValuev);
        _outFile << *it;
        while (++it != end(colValuev)) _outFile << ' ' << *it;
        _outFile << '\n';
    }
}
#endif  // CORE_SIMPLECOLUMNFILE_HPP
