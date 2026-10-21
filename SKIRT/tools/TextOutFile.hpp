/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TEXTOUTFILE_HPP
#define TEXTOUTFILE_HPP

#include "Array.hpp"
#include "CompileTimeUtils.hpp"
#include <array>
#include <fstream>
class FilePaths;
class Log;
class SimulationItem;
class Units;

////////////////////////////////////////////////////////////////////

/** This class allows writing text to a file specified in the constructor, with explicit support
    for formatting columns of floating point or integer numbers. Text is written per line, by
    calling the writeLine() or writeRow() functions. In a multiprocessing environment, only the
    root process will be allowed to write to the specified file; calls to writeLine() or writeRow()
    performed by other processes will have no effect.

    By default, the text is written to a regular file, regardless of the output path of the
    simulation. If the client specifies a text type in the constructor and the output path
    includes an HDF5 file, the text is instead collected in memory and written as an unstructured
    text file bundle in that HDF5 file when the file is closed. Such a bundle holds a scalar
    dataset \em text with the complete text as a string of variable length in UTF-8 encoding, and
    the attribute \em type with the specified text type, in addition to the standard attributes of
    every bundle (see the H5Lib class). The text types used by SKIRT are "free form" (for human
    readable text written by probes), "XML" (for the parameter file), and "log" (for the log
    file). In this case, the client must call the close() function explicitly, because writing a
    bundle may fail, and only that function can report the failure as a fatal error.

    The static copyToBundle() function stores a finished regular text file as an unstructured text
    file bundle, if the output path includes an HDF5 file. It is used for files that must remain
    available as regular files, such as the log file. */
class TextOutFile
{
    //=============== Construction - Destruction  ==================

public:
    /** The constructor of the TextOutFile class. If the constructor is invoked from the root
        process, the output file is opened and a stream for the file is initialized. For other
        processes, the output stream remains uninitialized and is never used, i.e. calling the
        writeLine() function has no effect. The constructor takes several arguments: (1) \em item
        specifies a simulation item in the hierarchy of the caller (usually the caller itself) used
        to retrieve the output file path and an appropriate logger, and to determine whether this
        is the root process; (2) \em filename specifies the name of the file, excluding path,
        simulation prefix and filename extension; (3) \em description describes the contents of the
        file for use in the log message issued after the file is successfully closed; and (4) the
        optional \em textType specifies the type of text for an unstructured text file bundle. If
        the text type is specified and the output path of the simulation includes an HDF5 file, no
        file is opened, and the text is collected for writing a bundle, as described in the class
        header. */
    TextOutFile(const SimulationItem* item, const string& filename, const string& description,
                const string& textType = string());

    /** In the root process, this function closes the file, or writes the bundle with the text
        collected so far, and logs an informational message, if this was not already done. Writing
        a bundle may fail, in which case the function throws a fatal error. It is important to call
        close() or allow the object to go out of scope before logging other messages or starting
        another significant chunk of work. */
    void close();

    /** The destructor does not throw exceptions. It calls the close() function, unless an
        exception is being propagated, in which case a bundle is not written. If close() throws a
        fatal error, the destructor logs the error message instead of propagating the error. */
    ~TextOutFile();

    /** If the output path of the simulation includes an HDF5 file, this function stores the
        contents of the regular text file with the specified name as an unstructured text file
        bundle in that HDF5 file, with the specified text type. Otherwise, the function does
        nothing. The arguments are (1) a simulation item in the hierarchy of the caller, used to
        retrieve the output file path and an appropriate logger; (2) the name of the file,
        including its filename extension but excluding path and simulation prefix, as for the
        FilePaths::output() function; (3) a description of the contents of the file, for use in
        error messages; and (4) the text type. The file must have been written completely, and it
        is left in place. The function returns the location of the bundle, or the empty string if
        it did not store a bundle. It does not log a message, so that it can store the log file of
        the simulation without changing it: tools that verify whether a simulation has finished
        expect a specific message at the end of that file. In a multiprocessing environment, only
        the root process stores the bundle. */
    static string copyToBundle(const SimulationItem* item, const string& name, const string& description,
                               const string& textType);

    //====================== Other functions =======================

public:
    /** This function writes the specified string to the file as a new line, or adds it to the text
        collected for a bundle. If the calling process is not the root, this function will have no
        effect. */
    void writeLine(const string& line);

    /** This function writes any buffered output to the regular file, so that the contents written
        so far can be inspected while the simulation progresses. It has no effect for a bundle, or
        if the calling process is not the root. */
    void flush();

    /** This function (virtually) adds a new column to the text file, characterized by a certain
        description and formatting. The format is 'd' for integer values, 'e' for scientific
        notation, 'f' for floating point notation and 'g' for the most concise 'e' or 'f'. For the
        'e' and 'f' formats, the precision represents the number of digits after the decimal point.
        For the 'g' format, the precision represents the maximum number of significant digits
        (trailing zeroes are omitted). The description of each column is added to the header of the
        text file, along with the column number. */
    void addColumn(const string& quantityDescription, string unitDescription = string(), char format = 'e',
                   int precision = 9);

    /** This function writes the specified list of (double) values to the text file, on a single row
        where adjacent values are seperated by a space. The values are formatted according to the
        'format' and 'precision' specified by the addColumn function. If the number of values in the
        list does not match the number of columns, a FatalError is thrown. */
    void writeRow(vector<double> values);

    /** This template function writes the specified list of values to the text file, on a single
        row where adjacent values are seperated by a space. The values are formatted according to
        the 'format' and 'precision' specified by the addColumn function. If the number of values
        in the list does not match the number of columns, a FatalError is thrown. */
    template<typename... Values, typename = std::enable_if_t<CompileTimeUtils::isNumericArgList<Values...>()>>
    void writeRow(Values... values)
    {
        std::array<double, sizeof...(values)> list = {{static_cast<double>(values)...}};
        writeRowPrivate(sizeof...(values), &list[0]);
    }

private:
    /** This function writes the specified list of (double) values to the text file with the same
        semantics as the other writeRow() functions. It is intended for private use from the
        template writeRow() functions. */
    void writeRowPrivate(size_t n, const double* values);

    //======================== Data Members ========================

protected:
    // can be used by subclasses
    Units* _units{nullptr};  // for conversion to output units
    std::ofstream _out;      // the output stream

private:
    // used for column formatting
    size_t _ncolumns{0};
    vector<char> _formats;
    vector<int> _precisions;

    // used when closing
    Log* _log{nullptr};  // the logger
    string _message;     // the message

    // used when writing a bundle
    FilePaths* _filePaths{nullptr};  // the file paths of the simulation, for verifying the I/O cutover
    string _bundlePath;              // the location of the bundle, or empty if a regular file is written
    string _textType;                // the type of text
    string _text;                    // the text collected so far
};

////////////////////////////////////////////////////////////////////

#endif
