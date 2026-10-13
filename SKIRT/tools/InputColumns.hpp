/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef INPUTCOLUMNS_HPP
#define INPUTCOLUMNS_HPP

#include "Array.hpp"
class Log;
class Units;

////////////////////////////////////////////////////////////////////

/** An instance of the InputColumns class keeps track of the columns of a table being read from an
    input file, independently of the format of that file. It is a helper for the TextInFile and
    ColumnInFile classes, which read the actual data from a column text file or from a text column
    file bundle in an HDF5 file, and it ensures that both formats are handled in the same way.

    The class distinguishes between the "physical" columns in the file, as described by the
    information in the file (the column information lines in the header of a text file, or the
    datasets in a bundle), and the "logical" columns expected by the program. The client first
    registers the physical columns, in the order in which they occur in the file, by calling the
    addFileColumn() function. It then calls the useColumns() function (optionally) and the
    addColumn() function for each logical column, in the same way as described for the functions
    with the same name in the TextInFile class. For each logical column, this class determines the
    corresponding physical column, verifies the units, and determines the conversion from the
    input units to internal SKIRT units. Finally, when reading the data, the client uses the
    logicalIndices() function to map physical to logical columns, and the converted() or
    convertColumn() functions to convert the values to internal units. */
class InputColumns
{
    //================= Construction - Destruction =================

public:
    /** The constructor receives the units system used to verify and convert the input units, and
        the logger used to log the information on each logical column, which may be the null
        pointer if this information should not be logged. The last argument describes the origin
        of the information on the physical columns, for use in error messages; for example "file
        header" or "bundle". */
    InputColumns(Units* units, Log* log, string source);

    //==================== Declaring the columns ===================

public:
    /** This function registers the next physical column in the file, with the specified name and
        unit string. The name is the description of the column given in the file, with consecutive
        white space characters replaced by a single space and without white space at the start
        and at the end. */
    void addFileColumn(const string& name, const string& unit);

    /** This function registers the next physical column in the file, with the specified name, for
        a file that does not specify the unit of the column. The values in the column are then
        assumed to be in the default unit specified by the program. */
    void addFileColumn(const string& name);

    /** This function specifies a mapping between the physical columns and the logical columns, as
        described for the TextInFile::useColumns() function. */
    void useColumns(string columns);

    /** This function adds the next logical column, as described for the TextInFile::addColumn()
        function. */
    void addColumn(const string& description, const string& quantity, const string& defaultUnit);

    //===================== Reading the columns ====================

public:
    /** This constant is used as an index value to indicate the absence of a column. */
    static constexpr size_t noIndex = 999999;

    /** This function returns the number of logical columns added so far. */
    size_t numLogicalColumns() const { return _numLogCols; }

    /** This function returns true if at least one of the logical columns is a virtual column
        containing zero values, i.e. a column that does not correspond to a physical column. */
    bool hasZeroColumns() const { return _haveZeroCols; }

    /** This function returns a list with an item for each physical column that must be read, in
        the order in which these columns occur in the file. The list is shorter than the number of
        physical columns if the last physical columns are not needed. Each item provides the
        zero-based index of the logical column that corresponds to the physical column, or the
        value InputColumns::noIndex if the physical column is not needed. */
    const vector<size_t>& logicalIndices() const { return _logColIndices; }

    /** This function converts the specified value, read from the physical column that corresponds
        to the logical column with the specified zero-based index, from input units to internal
        units, and returns the result. The last argument provides the values for the logical
        columns in the same row. It must already hold the converted value for each logical column
        that corresponds to a physical column preceding the one being converted, because the
        conversion of a "specific" quantity depends on the wavelength in the same row. */
    double converted(size_t logicalIndex, double value, const Array& row) const
    {
        const ColumnInfo& col = _colv[logicalIndex];
        if (col.convPower != 1.) value = pow(value, col.convPower);
        return value * (col.waveExponent ? pow(row[col.waveIndex], col.waveExponent) : col.convFactor);
    }

    /** This function converts all values in the specified array, read from the physical column
        that corresponds to the logical column with the specified zero-based index, from input
        units to internal units. The last argument provides the arrays with the values for each of
        the logical columns. It must already hold the converted values for each logical column
        that corresponds to a physical column preceding the one being converted. The function
        performs exactly the same calculation as the converted() function does for a single
        value. */
    void convertColumn(size_t logicalIndex, Array& values, const vector<Array>& columns) const;

    //======================== Private helpers ========================

private:
    /** This function returns the zero-based index in the column info list of the column that has
        a file name equal to the given name, or an error value if there is no such column or if
        there are multiple such columns. */
    size_t indexForName(const string& name) const;

    /** This function returns the logical index of the first logical column that is described as
        "wavelength" and that both logically and physically precedes the current column, or the
        error value if there is no such column. */
    size_t waveIndexForSpecificQuantity() const;

    //======================== Data Members ========================

private:
    Units* _units{nullptr};  // the units system
    Log* _log{nullptr};      // the logger for column information, or null if this information is not logged
    string _source;          // the origin of the information on the physical columns, for use in error messages

    // private type to store column info
    class ColumnInfo
    {
    public:
        size_t physColIndex{0};  // one-based physical index of this column in the file
        string title;            // description specified in the file, used to remap columns
        string description;      // official description provided by the program
        string quantity;         // quantity, provided by the program
        string unit;             // unit, provided by the program or specified in the file
        bool hasFileUnit{true};  // false if the file has information on this column but does not specify its unit
        double convFactor{1.};   // unit conversion factor from input to internal
        double convPower{1.};    // unit conversion power (exponent) from input to internal
        int waveExponent{0};     // wavelength exponent for converting "specific" quantities
        size_t waveIndex{0};     // zero-based logical index of wavelength column for converting "specific" quantities
    };

    bool _hasFileInfo{false};    // becomes true if the file has column info
    bool _hasProgInfo{false};    // becomes true if the program has added at least one column
    bool _doAutoCols{false};     // becomes true if the user requested automatic column assignment
    bool _allowZeroCols{false};  // becomes true if the automatic column assignment allows zero columns
    bool _haveZeroCols{false};   // becomes true if there is at least one (explicit or automatic) virtual zero column

    vector<ColumnInfo> _colv;  // info for each column, derived from file info and/or program info
    size_t _numLogCols{0};     // number of logical columns, or number of program columns added so far

    vector<size_t> _logColIndices;  // zero-based index into _colv for each physical column to be read
};

////////////////////////////////////////////////////////////////////

#endif
