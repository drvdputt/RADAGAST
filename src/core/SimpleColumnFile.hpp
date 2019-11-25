#ifndef CORE_SIMPLECOLUMNFILE_HPP
#define CORE_SIMPLECOLUMNFILE_HPP

#include <string>
#include <vector>

/** Class to read in data file with tab or space separated columns. TODO: merge with
    IOTools::ColumnFile to add write functionality, i.e. fill in and use commented out
    members. */
class SimpleColumnFile
{
public:
	/** Create a new instance, setting the file name. No file handle is opened yet */
	SimpleColumnFile(const std::string& fname) : _fname{fname} {}

	/** Overwrite the file on storage using the columns stored in this object */
	// void write() const;

	/** Overwrite the columns stored in this object by reading the file from storage. The
	    data from the first numCols columns will be stored. For speed, a guess for the number
	    of lines can be given. */
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

#endif // CORE_SIMPLECOLUMNFILE_HPP
