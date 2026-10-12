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
    functions throw a fatal error. As a result, the code in this build target and its clients can
    be compiled without any dependency on the HDF5 library.

    The functions in this namespace refer to an open HDF5 object (a group or a dataset) through an
    opaque 64-bit integer handle. A handle obtained from one of the open functions must eventually
    be released with the close() function. A handle remains valid until it is closed, even after
    the handle for the group from which it was obtained has been closed. The HDF5 library keeps a
    file open as long as there is an open handle for one of the objects in the file.

    Several functions take a \em context argument. This is a human-readable description of the
    object on which the function operates, used in the message of a fatal error.

    All functions in this namespace are thread-safe. The HDF5 library cannot be assumed to be
    thread-safe on its own, since that depends on how the installed library was built. All calls
    into the library are therefore protected by a single mutex. */
namespace H5Core
{
    /** This function returns true if the code was built with HDF5 support, and false otherwise. */
    bool available();

    /** This function opens the HDF5 file with the specified file path for reading, opens the group
        with the specified path inside that file, and returns a handle for the group. An empty
        group path refers to the root group of the file. The function throws a fatal error if the
        file or the group cannot be opened. */
    std::int64_t openGroup(const string& filePath, const string& groupPath);

    /** This function opens the group with the specified name inside the specified parent group
        and returns a handle for it. The function throws a fatal error if there is no such group. */
    std::int64_t openGroup(std::int64_t group, const string& name, const string& context);

    /** This function returns the names of the groups that are located directly inside the
        specified group, in alphabetical order. Other objects in the group, such as datasets, are
        ignored. */
    vector<string> groupNames(std::int64_t group, const string& context);

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
}

////////////////////////////////////////////////////////////////////

#endif
