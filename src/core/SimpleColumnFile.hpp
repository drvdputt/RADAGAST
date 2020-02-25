#ifndef CORE_SIMPLECOLUMNFILE_HPP
#define CORE_SIMPLECOLUMNFILE_HPP

#include "Error.hpp"
#include <fstream>
#include <string>
#include <vector>

/** Class to read in data file with tab or space separated columns. TODO: merge with OutColumnFile
    to add write functionality, i.e. fill in and use commented out members. */
class InColumnFile
{
public:
    /** Create a new instance, setting the file name. No file handle is opened yet */
    InColumnFile(const std::string& fname) : _fname{fname} {}

    /** Overwrite the file on storage using the columns stored in this object */
    // void write() const;

    /** Overwrite the columns stored in this object by reading the file from storage. The data from
        the first numCols columns will be stored. For speed, a guess for the number of lines can be
        given. */
    void read(int numCols, int reserveLines = 0);

    /** Get the data from one of the columns */
    const std::vector<double>& column(int i) { return _columnv[i]; }

    /** Remove any column stored in this object */
    // void clear() const;

    /** Add a column to the object */
    // void addColumn(const std::vector<double>& _data, const std::string& name = "");

private:
    std::string _fname;
    // std::vector<std::string> _namev;
    std::vector<std::vector<double>> _columnv;
};

class OutColumnFile
{
public:
    OutColumnFile(const std::string& filePath, const std::vector<std::string>& colNamev);
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

#endif  // CORE_SIMPLECOLUMNFILE_HPP
