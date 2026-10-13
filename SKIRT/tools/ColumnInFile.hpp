/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef COLUMNINFILE_HPP
#define COLUMNINFILE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
class H5BundleR;
class Log;
class SimulationItem;
class TextInFile;

////////////////////////////////////////////////////////////////////

/** This class allows reading a table of floating point values from a user input file, with support
    for unit conversions. The table can be provided as a regular column text file, or as a text
    column file bundle in an HDF5 file. The class offers the same functions as the TextInFile
    class, so that a client reads both formats in the same way and does not need to know which of
    them is being used.

    The filename is specified in the constructor, and information about the columns expected in the
    file is specified by calling the addColumn() function for each column, possibly preceded by a
    call to the useColumns() function. Nothing is opened until the first row is read. At that
    time, the filename is passed to the FilePaths::input() function to obtain the candidate
    locations of the input file, taking into account the input path specified on the command line:

    - If there is a regular file at the candidate location for a regular file, that file is read
      through an internally held TextInFile object. See the TextInFile class for a description of
      the format of a column text file.

    - Otherwise, if there is a candidate location for a bundle (because the input path or the
      filename itself refers to an HDF5 file), that bundle is read as described below.

    If the input file is found at neither of two candidate locations, the fatal error that reports
    this mentions both the regular file and the bundle.

    This class is intended for user input files only; a resource file is never provided as a bundle
    and must be read through the TextInFile class.

    Previous versions of SKIRT also accepted a table in the binary "stored columns" format,
    recognized by the ".scol" filename extension. This format is no longer supported. A filename
    with that extension, whether it refers to a regular file or to a bundle, causes a fatal error
    when the first row is read. Such a file must be converted to a regular column text file or to
    a text column file bundle, for example with the conversion commands offered by PTS.

    Text column file bundle
    -----------------------

    A text column file bundle holds a one-dimensional dataset for each column of the table. All
    datasets must have the same length, the number of rows. A dataset is named for the column it
    represents, i.e. it has the description that would be listed for the column in the header of a
    column text file, and it has the following attributes:

    - \em column: the zero-based index of the column in the table (an integer). The values of this
      attribute for the datasets in a bundle must be consecutive and start at zero. They define the
      order of the columns, which is otherwise undefined because the datasets in a bundle are not
      ordered. Note that the column numbers in the header of a column text file are one-based.

    - \em unit: the unit string for the values in the column, exactly as in the header of a column
      text file, with "1" or the empty string for a dimensionless quantity. If this attribute is
      missing, the values are assumed to be in the default units specified by the program, as for
      a column text file without column information.

    Datasets without a \em column attribute are ignored. Because a forward slash cannot occur in
    the name of a dataset, a vertical bar in the name of a dataset stands for a forward slash in
    the description of the column; for example, the column "n_He/n_H" is stored as the dataset
    "n_He|n_H".

    The columns in a bundle are matched to the columns expected by the program in exactly the same
    way as the columns in a column text file: by default in the order given by the \em column
    attribute, or by name as requested through the useColumns() function. Only the datasets for the
    columns actually needed are read, each one as a whole. The values are converted to the 64-bit
    floating point type from the type in which they are stored, and the data may be chunked or
    compressed. */
class ColumnInFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor records its arguments without opening anything. The arguments are the same
        as for the TextInFile constructor: (1) \em item specifies a simulation item in the
        hierarchy of the caller (usually the caller itself) used to retrieve the input file path,
        the units system, and an appropriate logger; (2) \em filename specifies the name of the
        file as provided by the user, including filename extension but excluding path; (3) \em
        description describes the contents of the file for use in log messages. */
    ColumnInFile(const SimulationItem* item, const string& filename, const string& description);

    /** The copy constructor is deleted because instances of this class cannot be copied. */
    ColumnInFile(const ColumnInFile&) = delete;

    /** The assignment operator is deleted because instances of this class cannot be copied. */
    ColumnInFile& operator=(const ColumnInFile&) = delete;

    /** This function closes the file or releases the data read from the bundle, if this was not
        already done. It is best to call close() or allow the object to go out of scope before
        logging other messages or starting another significant chunk of work. */
    void close();

    /** The destructor calls the close() function. */
    virtual ~ColumnInFile();

    //====================== Declaring columns =======================

public:
    /** This function specifies a mapping between the "physical" columns in the file and the
        "logical" columns expected by the program, as described for the TextInFile::useColumns()
        function. It must be called before the first invocation of the addColumn() function. For
        a bundle, the name of a physical column is the name of the corresponding dataset, with any
        vertical bar replaced by a forward slash. */
    void useColumns(const string& columns);

    /** This function adds a column to the list of columns expected by the program, characterized
        by the given description, physical quantity, and default unit, as described for the
        TextInFile::addColumn() function. The default unit is used for a column text file without
        column information, and for a dataset in a bundle that does not specify its unit. All
        columns must be added before the first row is read. */
    void addColumn(const string& description, const string& quantity = string(), const string& defaultUnit = string());

    //========================== Reading ===========================

public:
    /** This function reads the next row and stores the resulting values, converted to internal
        units, in the array passed to the function by reference. The array is resized to the
        number of columns added by the program. The function returns true if a row was read, and
        false if there are no more rows. If no columns have been added, the function skips the
        next row and sets the array to zero length. See the TextInFile::readRow() function for
        more information. */
    bool readRow(Array& values);

    /** This variadic template function reads the next row and stores the resulting values in the
        variables passed to the function by reference. The number of arguments must match the
        number of columns added by the program. For example:

        \verbatim
        double a,b,c,d;
        bool success = in.readRow(a,b,c,d);  // reads a row from a file with 4 columns
        \endverbatim
        */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isFloatArgList<Values...>()>>
    bool readRow(Values&... values)
    {
        Array result;
        bool success = readRow(result);
        if (success) assignValues(0, result, values...);
        return success;
    }

    /** This function reads all remaining rows and returns the resulting values as a vector of row
        arrays. For each row, this function behaves just like readRow(Array&). */
    vector<Array> readAllRows();

    /** This function reads all remaining rows and returns the resulting values as a vector of
        column arrays, with a column for each of the columns added by the program. */
    vector<Array> readAllColumns();

    /** This function reads all remaining rows and stores the resulting column arrays in the
        variables passed to the function by reference. The number of arguments must match the
        number of columns added by the program. */
    template<typename... Columns, typename = std::enable_if_t<sizeof...(Columns) != 0>>
    void readAllColumns(Columns&... columns)
    {
        auto result = readAllColumns();
        assignColumns(0, result, columns...);
    }

    //================= Functions for use in subclasses =================

protected:
    /** This function locates and opens the input file, as described in the header of this class,
        and processes the columns that were declared, unless this already happened. It is called
        by all reading functions in this class, and must be called by a subclass before it reads
        information on its own. */
    void open();

    /** This function returns the internally held TextInFile object if the input is being read from
        a regular file, and the null pointer if it is being read from a bundle. It must be called
        after the open() function. */
    TextInFile* textFile() const { return _text.get(); }

    /** This function is called by the open() function when the input is being read from a bundle,
        just after the bundle has been opened and before any of its column datasets are read. The
        default implementation does nothing. A subclass that supports a bundle with extra datasets
        can override this function to read these datasets, and to call the setNumRows() function
        if the extra datasets determine the number of rows. */
    virtual void openedBundle(const H5BundleR& bundle);

    /** This function specifies the number of rows in the bundle being read. The open() function
        then verifies that each of the column datasets has this number of rows, and reports this
        number of rows if no columns are being read. It may be called only from the
        openedBundle() function. */
    void setNumRows(size_t numRows);

    /** This function returns the number of rows read from a bundle so far. */
    size_t numRowsRead() const { return _nextRow; }

    /** This function returns the location of the bundle being read, for use in error messages. It
        must be called after the open() function. */
    const string& bundlePath() const { return _bundlePath; }

    //======================== Private helpers ========================

private:
    // opens the bundle with the specified location and loads the columns that were declared
    void openBundle();

    // recursively assign values from Array to double& arguments; used in variadic readRow()
    template<typename... Values>
    static inline void assignValues(size_t index, const Array& result, double& value, Values&... values)
    {
        value = result[index];
        assignValues(index + 1, result, values...);
    }
    static inline void assignValues(size_t /*index*/, const Array& /*result*/) {}

    // recursively assign columns from vector to Array& arguments; used in variadic readAllColumns()
    template<typename... Columns>
    static inline void assignColumns(size_t index, vector<Array>& result, Array& column, Columns&... columns)
    {
        column = std::move(result[index]);
        assignColumns(index + 1, result, columns...);
    }
    static inline void assignColumns(size_t /*index*/, vector<Array>& /*result*/) {}

    //======================== Data Members ========================

private:
    // the information recorded by the constructor and by the functions that declare columns
    const SimulationItem* _item{nullptr};
    string _filename;
    string _description;
    struct Declaration
    {
        bool isUseColumns{false};      // true for a call to useColumns(), false for a call to addColumn()
        string first, quantity, unit;  // the arguments of the call; "first" holds the columns or the description
    };
    vector<Declaration> _declarations;

    // the state of the file
    bool _isOpen{false};    // becomes true when the file or the bundle has been opened
    bool _isClosed{false};  // becomes true when the file or the bundle has been closed

    // the regular file, if the input is being read from a regular file
    std::unique_ptr<TextInFile> _text;

    // the data read from the bundle, if the input is being read from a bundle
    Log* _log{nullptr};       // the logger
    string _bundlePath;       // the location of the bundle
    vector<Array> _columns;   // the values for each logical column, in internal units
    size_t _numLogCols{0};    // the number of logical columns
    size_t _numRows{0};       // the number of rows
    size_t _nextRow{0};       // the zero-based index of the next row to be read
    bool _hasNumRows{false};  // becomes true when the number of rows is known
};

////////////////////////////////////////////////////////////////////

#endif
