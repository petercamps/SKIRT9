/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "TextInFile.hpp"
#include "FatalError.hpp"
#include "FilePaths.hpp"
#include "Log.hpp"
#include "StringUtils.hpp"
#include "System.hpp"
#include "Units.hpp"
#include <exception>
#include <regex>
#include <sstream>

////////////////////////////////////////////////////////////////////

namespace
{
    // This function looks for the next header line that conforms to the required structured syntax. If such a line
    // is found, the function stores the column index (or if missing, -1), description and unit string in the arguments
    // and returns true. If no such header line is found, the function consumes the complete header and returns false.
    bool getNextInfoLine(std::ifstream& in, int& colIndex, string& description, string& unit)
    {
        // continue reading until a conforming header line is found or until the complete header has been consumed
        while (true)
        {
            // consume whitespace characters but nothing else
            while (true)
            {
                auto ch = in.peek();
                if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r') break;
                in.get();
            }

            // if the first non-whitespace character is not a hash character, there is no header line
            if (in.peek() != '#') return false;

            // read the header line
            string line;
            getline(in, line);

            // if the line conforms to the required syntax, return the extracted information
            static const std::regex syntax("#\\s*column\\s*(\\d*)\\s*:\\s*([^()]*)\\(\\s*([a-zA-Z0-9/]*)\\s*\\)\\s*",
                                           std::regex::icase);
            std::smatch matches;
            if (std::regex_match(line, matches, syntax) && matches.size() == 4)
            {
                string index = matches[1].str();
                colIndex = index.empty() ? -1 : std::stoul(index);
                description = StringUtils::squeeze(matches[2].str());
                unit = matches[3].str();
                return true;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////

TextInFile::TextInFile(const SimulationItem* item, const string& filename, const string& description, bool resource,
                       bool silent)
    : _log(item->find<Log>()), _columns(item->find<Units>(), resource ? nullptr : _log, "file header"),
      _isResource(resource)
{
    // get the full path for the resource or for the input file; an input file name that refers to a bundle in an
    // HDF5 file has no regular file candidate
    string filepath = resource ? FilePaths::resource(filename) : item->find<FilePaths>()->input(filename).first;
    if (filepath.empty())
        throw FATALERROR("A column text file cannot be read from a bundle in an HDF5 file: " + filename);

    // open the file
    _in = System::ifstream(filepath);
    if (!_in) throw FATALERROR("Could not open the " + description + " text file " + filepath);

    // log "reading file" message
    if (!silent) _log->info(item->typeAndName() + " reads " + description + " from text file " + filepath + "...");

    // register the descriptions and units of the columns specified in any structured header lines
    int index;  // one-based column index obtained from file info
    string title;
    string unit;
    size_t numFileColumns = 0;
    while (getNextInfoLine(_in, index, title, unit))
    {
        numFileColumns++;
        if (index != -1 && static_cast<size_t>(index) != numFileColumns)
            throw FATALERROR("Incorrect column index in file header for column " + std::to_string(numFileColumns));
        _columns.addFileColumn(title, unit);
    }
    _hasTextOpen = true;
}

////////////////////////////////////////////////////////////////////

void TextInFile::close()
{
    if (_hasTextOpen)
    {
        _in.close();
        _hasTextOpen = false;

        // log "done" message, except if this is a resource file or after an exception has been thrown
        if (!_isResource && std::uncaught_exceptions() == 0) _log->info("Done reading");
    }
}

////////////////////////////////////////////////////////////////////

TextInFile::~TextInFile()
{
    close();
}

////////////////////////////////////////////////////////////////////

void TextInFile::useColumns(string columns)
{
    _columns.useColumns(columns);
}

////////////////////////////////////////////////////////////////////

void TextInFile::addColumn(const string& description, const string& quantity, string defaultUnit)
{
    _columns.addColumn(description, quantity, defaultUnit);
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readRow(Array& values)
{
    // read next row
    if (_hasTextOpen)
    {
        // read new line until it is non-empty and non-comment
        string line;
        while (_in.good())
        {
            getline(_in, line);
            if (auto pos = line.find_first_not_of(" \t"); pos != string::npos && line[pos] != '#')
            {
                // resize result array (and clear it in case there are virtual zero columns)
                size_t numLogCols = _columns.numLogicalColumns();
                if (values.size() != numLogCols || _columns.hasZeroColumns()) values.resize(numLogCols);

                // convert values from line and store them in result array
                std::stringstream linestream(line);
                for (size_t i : _columns.logicalIndices())  // i: zero-based logical index
                {
                    if (linestream.eof()) throw FATALERROR("One or more required value(s) on text line are missing");

                    // read the value as floating point
                    double value;
                    linestream >> value;
                    if (linestream.fail())
                    {
                        // some compilers/libraries do not support reading NaN or Inf values, while others do;
                        // we here provide backup support for NaN values (infinities are more complex because of
                        // the need for handling the negative sign which may be already consumed by the stream)
                        linestream.clear();
                        string offending;
                        linestream >> offending;
                        if (StringUtils::toLower(offending) != "nan")
                            throw FATALERROR("Input text is not formatted as a floating point number: " + offending);
                        value = std::numeric_limits<double>::quiet_NaN();
                    }

                    // if mapped to a logical column, convert from input units to internal units, and store the result
                    if (i != InputColumns::noIndex) values[i] = _columns.converted(i, value, values);
                }
                return true;
            }
        }
    }

    // end of file was reached or no file is open
    return false;
}

////////////////////////////////////////////////////////////////////

bool TextInFile::readNonLeaf(int& nx, int& ny, int& nz)
{
    // read next non-leaf row
    if (_hasTextOpen)
    {
        string line;

        while (true)
        {
            int c = _in.peek();

            // skip comments line
            if (c == '#')
            {
                getline(_in, line);
            }

            // process nonleaf line
            else if (c == '!')
            {
                _in.get();  // skip exclamation mark
                getline(_in, line);

                // convert nx,ny,nz values from line and store them in output arguments
                std::stringstream linestream(line);
                linestream >> nx >> ny >> nz;
                if (linestream.fail())
                    throw FATALERROR("Nonleaf subdivision specifiers are missing or not formatted as integers");

                return true;
            }

            // eat leading white space and empty lines
            else if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
            {
                _in.get();
            }

            // signal not a nonleaf line
            else
            {
                return false;
            }
        }
    }

    // no file is open
    return false;
}

////////////////////////////////////////////////////////////////////

vector<Array> TextInFile::readAllRows()
{
    vector<Array> rows;
    while (true)
    {
        rows.emplace_back();        // add a default-constructed array to the vector
        if (!readRow(rows.back()))  // read next line's values into that array
        {
            rows.pop_back();  // at the end, remove the extraneous array
            break;
        }
    }
    return rows;
}

////////////////////////////////////////////////////////////////////

vector<Array> TextInFile::readAllColumns()
{
    // read the remainder of the file into rows
    const vector<Array>& rows = readAllRows();
    size_t nrows = rows.size();
    size_t ncols = _columns.numLogicalColumns();

    // transpose the result into columns
    vector<Array> columns(ncols, Array(nrows));
    for (size_t c = 0; c != ncols; ++c)
        for (size_t r = 0; r != nrows; ++r) columns[c][r] = rows[r][c];

    return columns;
}

////////////////////////////////////////////////////////////////////
