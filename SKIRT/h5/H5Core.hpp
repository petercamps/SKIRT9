/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5CORE_HPP
#define H5CORE_HPP

#include "Array.hpp"

////////////////////////////////////////////////////////////////////

/** The H5Core namespace is the interface to the HDF5 C library, and its implementation is the only
    place in the SKIRT code where explicit calls to this library are allowed. It is intended for
    use by the classes in the \c h5 build target only; all other code accesses HDF5 files
    through the H5Lib class and the classes it hands out.

    If the code is built without HDF5 support (i.e. the \c BUILD_HDF5 build option is disabled),
    the available() function returns false, the close() function does nothing, and all other
    functions, for reading and for writing, throw a fatal error. As a result, the code in this
    build target and its clients can be compiled without any dependency on the HDF5 library.

    The functions in this namespace refer to an open HDF5 object (a group or a dataset) through an
    opaque 64-bit integer handle. A handle obtained from one of the open functions must eventually
    be released with the close() function. A handle remains valid until it is closed, even after
    the handle for the group from which it was obtained has been closed. The HDF5 library keeps a
    file open as long as there is an open handle for one of the objects in the file.

    Several functions take a \em context argument. This is a human-readable description of the
    object on which the function operates, used in the message of a fatal error.

    The functions for writing store each value in the data type listed for it, without any
    conversion: a Boolean as an enumeration based on a signed byte with the members \c FALSE (0)
    and \c TRUE (1), the representation used by the \c h5py Python package; an integer as a 32-bit
    signed integer; a floating point number as a 64-bit floating point number; and a string as a
    string of variable length in UTF-8 encoding. Datasets are stored contiguously, without
    chunking or compression.

    All functions in this namespace are thread-safe. The HDF5 library cannot be assumed to be
    thread-safe on its own, since that depends on how the installed library was built. All calls
    into the library are therefore protected by a single mutex.

    Within a process, the HDF5 library shares a single object for each open file. It cannot reopen
    a file for writing while the file is open for reading only, but it can reopen a file for
    reading while the file is open for writing. When several simulations run in parallel in the
    same process, one of them may hold a file open for reading while it loads its input, while
    another one is about to write its output to the same file. The functions that open a file for
    writing, createGroup() and verifyWritable(), therefore wait as long as any object in this
    process holds the file open for reading only. Readers hold a file only while they load data, so
    that the wait is short. If the file remains held for more than ten minutes, for example because
    the calling thread holds it itself, the function throws a fatal error rather than waiting
    forever. */
namespace H5Core
{
    /** This function returns true if the code was built with HDF5 support, and false otherwise. */
    bool available();

    /** This function opens the HDF5 file with the specified file path for reading, opens the group
        with the specified path inside that file, and returns a handle for the group. An empty
        group path refers to the root group of the file. The function throws a fatal error if the
        file or the group cannot be opened. */
    std::int64_t openGroup(const string& filePath, const string& groupPath);

    /** This function returns true if the HDF5 file with the specified file path exists and holds
        a group with the specified path, and false if the file does not exist or does not hold
        such a group. The function throws a fatal error if the file exists but cannot be opened as
        an HDF5 file. */
    bool hasGroup(const string& filePath, const string& groupPath);

    /** This function opens the group with the specified name inside the specified parent group
        and returns a handle for it. The function throws a fatal error if there is no such group. */
    std::int64_t openGroup(std::int64_t group, const string& name, const string& context);

    /** This function returns the names of the groups that are located directly inside the
        specified group, in alphabetical order. Other objects in the group, such as datasets, are
        ignored. */
    vector<string> groupNames(std::int64_t group, const string& context);

    /** This function returns the names of the datasets that are located directly inside the
        specified group, in alphabetical order. Other objects in the group, such as groups, are
        ignored. */
    vector<string> datasetNames(std::int64_t group, const string& context);

    /** This function returns true if the specified group directly contains a dataset with the
        specified name, and false otherwise. */
    bool hasDataset(std::int64_t group, const string& name);

    /** This function opens the dataset with the specified name inside the specified group and
        returns a handle for it. The function throws a fatal error if there is no such dataset. */
    std::int64_t openDataset(std::int64_t group, const string& name, const string& context);

    /** This function releases the specified handle. It does nothing if the handle is negative, and
        it never throws an exception. */
    void close(std::int64_t object) noexcept;

    /** This function returns true if the specified object (a group or a dataset) has an attribute
        with the specified name, and false otherwise. */
    bool hasAttribute(std::int64_t object, const string& name);

    /** This function returns the value of the Boolean attribute with the specified name of the
        specified object. The attribute must hold a single value of an integer or enumeration
        type; a nonzero value means true. Otherwise the function throws a fatal error. */
    bool boolAttribute(std::int64_t object, const string& name, const string& context);

    /** This function returns the value of the integer attribute with the specified name of the
        specified object. The attribute must hold a single value of an integer or floating point
        type that can be represented as a 32-bit signed integer without loss of information.
        Otherwise the function throws a fatal error. */
    int intAttribute(std::int64_t object, const string& name, const string& context);

    /** This function returns the value of the floating point attribute with the specified name of
        the specified object. The attribute must hold a single value of an integer or floating
        point type. Otherwise the function throws a fatal error. */
    double doubleAttribute(std::int64_t object, const string& name, const string& context);

    /** This function returns the value of the string attribute with the specified name of the
        specified object, as a UTF-8 encoded string. The attribute must hold a single string of
        fixed or variable length. Otherwise the function throws a fatal error. */
    string stringAttribute(std::int64_t object, const string& name, const string& context);

    /** This function returns the shape of the specified dataset, i.e. the number of elements along
        each of its dimensions, with the dimension that varies least rapidly in the stored order
        listed first. The result is empty for a dataset holding a single scalar value. */
    vector<size_t> shape(std::int64_t dataset, const string& context);

    /** This function returns the complete contents of the specified dataset as Boolean values, in
        the stored order, i.e. with the last dimension varying most rapidly. Each element of the
        result is 1 for true or 0 for false. The dataset must have an integer or enumeration type;
        a nonzero stored value means true. Otherwise the function throws a fatal error. */
    vector<char> readBool(std::int64_t dataset, const string& context);

    /** This function returns the complete contents of the specified dataset as 32-bit signed
        integers, in the stored order, i.e. with the last dimension varying most rapidly. The
        dataset must have an integer or floating point type, and each of its values must be
        representable without loss of information. Otherwise the function throws a fatal error. */
    vector<int> readInt(std::int64_t dataset, const string& context);

    /** This function returns the complete contents of the specified dataset as 64-bit floating
        point numbers, in the stored order, i.e. with the last dimension varying most rapidly. The
        dataset must have an integer or floating point type, and each of its values must be
        representable without loss of information. Otherwise the function throws a fatal error. */
    Array readDouble(std::int64_t dataset, const string& context);

    /** This function opens the HDF5 file with the specified file path for writing, creating the
        file if it does not exist, and creates a new group with the specified path inside that
        file, including any missing parent groups. It returns a handle for the new group. If the
        file already holds a group with that path that does not contain other groups (a bundle),
        that group is removed first; the rest of the file remains untouched. The function throws a
        fatal error if the file cannot be opened or created for writing, if the path refers to a
        dataset or to a group that contains other groups (a suite), or if one of the parent groups
        is a dataset or a group that directly contains datasets (a bundle). */
    std::int64_t createGroup(const string& filePath, const string& groupPath);

    /** This function sets the Boolean attribute with the specified name of the specified object to
        the specified value, replacing an existing attribute with the same name. */
    void setBoolAttribute(std::int64_t object, const string& name, bool value, const string& context);

    /** This function sets the integer attribute with the specified name of the specified object to
        the specified value, replacing an existing attribute with the same name. */
    void setIntAttribute(std::int64_t object, const string& name, int value, const string& context);

    /** This function sets the floating point attribute with the specified name of the specified
        object to the specified value, replacing an existing attribute with the same name. */
    void setDoubleAttribute(std::int64_t object, const string& name, double value, const string& context);

    /** This function sets the string attribute with the specified name of the specified object to
        the specified UTF-8 encoded value, replacing an existing attribute with the same name. */
    void setStringAttribute(std::int64_t object, const string& name, const string& value, const string& context);

    /** This function creates a dataset with the specified name and shape inside the specified
        group, writes the specified Boolean values to it, and returns a handle for the dataset.
        Each value is a character holding 1 for true or 0 for false. The values are listed with the
        last dimension varying most rapidly, and their number must equal the product of the
        dimensions in the shape. The function throws a fatal error if the dataset cannot be created
        or written, for example because the group already holds an object with that name. */
    std::int64_t createBoolDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                   const char* values, const string& context);

    /** This function creates a dataset with the specified name and shape inside the specified
        group, writes the specified 32-bit integer values to it, and returns a handle for the
        dataset, as described for the createBoolDataset() function. */
    std::int64_t createIntDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                  const int* values, const string& context);

    /** This function creates a dataset with the specified name and shape inside the specified
        group, writes the specified 64-bit floating point values to it, and returns a handle for
        the dataset, as described for the createBoolDataset() function. */
    std::int64_t createDoubleDataset(std::int64_t group, const string& name, const vector<size_t>& shape,
                                     const double* values, const string& context);

    /** This function creates a dataset with the specified name inside the specified group, holding
        the specified UTF-8 encoded text as a single string of variable length, and returns a
        handle for the dataset. The function throws a fatal error if the dataset cannot be created
        or written. */
    std::int64_t createTextDataset(std::int64_t group, const string& name, const string& text, const string& context);

    /** This function opens the HDF5 file with the specified file path for writing, creating the
        file if it does not exist, and closes it again. It throws a fatal error if the file cannot
        be opened or created for writing, as for the createGroup() function. It allows a client to
        verify that a file can be written before it writes to the file. */
    void verifyWritable(const string& filePath);

    /** This function returns the paths of the HDF5 files on which this process holds open
        objects, in the form passed to the function that opened them, without duplicates. If the
        code was built without HDF5 support, the list is empty. */
    vector<string> openFiles();
}

////////////////////////////////////////////////////////////////////

#endif
