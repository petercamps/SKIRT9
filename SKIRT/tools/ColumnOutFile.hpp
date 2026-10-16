/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COLUMNOUTFILE_HPP
#define COLUMNOUTFILE_HPP

#include "CompileTimeUtils.hpp"
#include <array>
class SimulationItem;
class TextOutFile;

////////////////////////////////////////////////////////////////////

/** This class allows writing a table of numbers to an output file, with a row for each item (for
    example a wavelength, a spatial cell, or a position) and a column for each property. The table
    is written as a regular column text file, or as a text column file bundle in an HDF5 file,
    depending on the output path of the simulation. The class offers the same functions for both
    formats, so that a client writes the table in the same way and does not need to know which of
    them is being used.

    The client specifies the file name and a description in the constructor, optionally a long
    description through the setLongDescription() function, and the columns through the
    addColumn() function. It then writes the rows through the writeRow() function, and finally
    calls the close() function. All columns must be added before the first row is written.
    Nothing is opened until the first row is written, or until the file is closed if there are no
    rows. At that time, the filename is passed to the FilePaths::output() function to obtain the
    candidate locations of the output file:

    - If the output path of the simulation includes an HDF5 file, the table is written as a bundle
      in that file, as described below. The rows are collected in memory, and the bundle is written
      when the file is closed.

    - Otherwise, the table is written as a regular column text file through an internally held
      TextOutFile object, row by row, exactly as a client that uses the TextOutFile class directly
      would write it: the long description on the first line, preceded by a hash character; a line
      with column information for each column, of the form "# column N: description (unit)"; and
      a line for each row, with the values formatted as specified for each column.

    A client must call the close() function explicitly, because writing a bundle may fail, and
    only that function can report the failure as a fatal error. The destructor completes a file
    that has not been closed, but because a destructor must not throw an exception, it can only
    log an error message, after which the simulation continues.

    In a multiprocessing environment, only the root process writes the output file; the functions
    called by other processes have no effect, except for reporting errors in their arguments.

    Text column file bundle
    -----------------------

    A text column file bundle holds a one-dimensional dataset of 64-bit floating point numbers for
    each column, with the number of rows as its length. The dataset has the attribute \em column,
    the zero-based index of the column, and the attribute \em unit, the unit string of the column,
    with "1" for a dimensionless quantity. The format and precision specified for a column apply to
    a text file only, because a bundle stores the values in binary form.

    The name of a dataset is derived from the description of its column by the following rules,
    applied in this order: white space at the start and the end is removed, and consecutive white
    space characters are replaced by a single space; a forward slash is replaced by a vertical
    bar, because HDF5 uses the forward slash to separate the names of nested groups; an empty
    description is replaced by the word "column" followed by the zero-based index of the column;
    and if the resulting name is already used for a preceding column, a sequence number between
    parentheses is added, as in "weight (2)". These are the rules of the data model for text
    column file bundles, which the \c pts.hdf5 package in PTS implements as well.

    The bundle has the attribute \em description, holding the long description, followed on
    separate lines by any comments added through the addComment() function, if there is any such
    text, in addition to the standard attributes of every bundle (see the H5Lib class).

    Incrementally written files
    ---------------------------

    A file that must be written row by row while the simulation progresses, so that its contents
    can be inspected during the run, is written through the TextOutFile class as a regular file,
    regardless of the output path. Once the file is complete, the static copyToBundle() function
    stores its contents in a text column file bundle, if the output path includes an HDF5 file. */
class ColumnOutFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor records its arguments without opening anything: (1) \em item specifies a
        simulation item in the hierarchy of the caller (usually the caller itself) used to retrieve
        the output file path and an appropriate logger; (2) \em filename specifies the name of the
        file, excluding path, simulation prefix and filename extension; (3) \em description
        describes the contents of the file for use in the log message issued after the file has
        been written. */
    ColumnOutFile(const SimulationItem* item, const string& filename, const string& description);

    /** The copy constructor is deleted because instances of this class cannot be copied. */
    ColumnOutFile(const ColumnOutFile&) = delete;

    /** The assignment operator is deleted because instances of this class cannot be copied. */
    ColumnOutFile& operator=(const ColumnOutFile&) = delete;

    /** This function completes the output file, if this was not already done. For a regular file,
        it closes the file. For a bundle, it writes the bundle with all rows written so far, and
        throws a fatal error if this fails. In both cases, it logs an informational message. A
        client must call this function after writing the last row, before logging other messages
        or starting another significant chunk of work. */
    void close();

    /** The destructor does not throw exceptions. If the client did not call the close() function,
        the destructor calls it, unless an exception is being propagated, in which case the output
        file is left incomplete. If close() throws a fatal error, the destructor logs the error
        message instead of propagating the error. */
    ~ColumnOutFile();

    //====================== Declaring the header ======================

public:
    /** This function specifies a free-form description of the contents of the file, which is
        written on the first line of a regular file, preceded by a hash character, and as the
        description attribute of a bundle. It must be called before the first row is written. */
    void setLongDescription(const string& longDescription);

    /** This function adds a column to the table, characterized by the specified description, unit
        string, and formatting, as described for the TextOutFile::addColumn() function. The format
        and precision apply only to a regular file. Columns must be added before the first row is
        written. */
    void addColumn(const string& quantityDescription, string unitDescription = string(), char format = 'e',
                   int precision = 9);

    /** This function adds a comment that further describes the contents of the file. In a regular
        file, the comment is written as a line in the header, preceded by a hash character, after
        the column information lines for the columns that have been added before the call. In a
        bundle, the comment becomes a separate line of the description attribute, after the long
        description. It must be called before the first row is written. */
    void addComment(const string& comment);

    //========================== Writing ===========================

public:
    /** This function writes a row with the specified values, one for each column, in the order in
        which the columns have been added. If the number of values does not match the number of
        columns, a fatal error is thrown. */
    void writeRow(const vector<double>& values);

    /** This template function writes a row with the specified values, one for each column, with
        the same semantics as the other writeRow() function. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isNumericArgList<Values...>()>>
    void writeRow(Values... values)
    {
        std::array<double, sizeof...(values)> list = {{static_cast<double>(values)...}};
        writeRowPrivate(sizeof...(values), &list[0]);
    }

    //=================== Incrementally written files ===================

public:
    /** If the output path of the simulation includes an HDF5 file, this function stores the
        contents of the regular column text file with the specified name, which has been written
        completely through the TextOutFile class, as a text column file bundle in that HDF5 file.
        Otherwise, the function does nothing. The arguments are the same as for the constructor of
        this class.

        The function parses the header of the file for the column descriptions and units, as the
        TextInFile class does for an input file, and creates the bundle and its datasets as an
        instance of this class would have created them. The comment lines in the header, other than
        the column information lines, become the description attribute of the bundle. The values
        on each data line are read as floating point numbers, one for each column. The regular file
        is left in place. In a multiprocessing environment, only the root process stores the
        bundle. */
    static void copyToBundle(const SimulationItem* item, const string& filename, const string& description);

    //======================== Private helpers ========================

private:
    /** This function determines the format of the output file and prepares it for writing, if
        this did not already happen. */
    void open();

    /** This function writes a row with the specified number of values, with the semantics of the
        public writeRow() functions. */
    void writeRowPrivate(size_t n, const double* values);

    /** This function writes a text column file bundle at the specified location, with the
        specified description attribute, and with a dataset for each column, named for the
        specified column description and holding the specified values and unit. It logs an
        informational message on behalf of the specified simulation item, mentioning the specified
        file description. The column values are released while writing the bundle. */
    static void writeBundle(const SimulationItem* item, const string& bundlePath, const string& description,
                            const string& bundleDescription, const vector<string>& descriptions,
                            const vector<string>& units, vector<vector<double>>& columns);

    //======================== Data Members ========================

private:
    // the information recorded by the constructor and by the functions that declare the header
    const SimulationItem* _item{nullptr};
    string _filename;
    string _description;
    string _longDescription;
    struct HeaderLine
    {
        bool isColumn{false};  // true for a column, false for a comment
        string text;           // the description of the column, or the comment
        string unit;           // the unit of the column
        char format{'e'};      // the format of the column
        int precision{9};      // the precision of the column
    };
    vector<HeaderLine> _header;
    size_t _numColumns{0};

    // the state of the file
    bool _isOpen{false};    // becomes true when the format of the output file has been determined
    bool _isClosed{false};  // becomes true when the output file has been completed

    // the regular file, if the table is being written as a regular file
    std::unique_ptr<TextOutFile> _text;

    // the location of the bundle, and the values collected for each column, if the table is being written as a bundle
    string _bundlePath;
    vector<vector<double>> _columns;
};

////////////////////////////////////////////////////////////////////

#endif
