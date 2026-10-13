/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef H5LIB_HPP
#define H5LIB_HPP

#include "H5SuiteR.hpp"

////////////////////////////////////////////////////////////////////

/** The H5Lib class offers the entry points for accessing HDF5 files through a thin C++ wrapper
    around the HDF5 C library, offering just the features used by the code. All code outside of
    the \c h5 build target depends only on this wrapper, never directly on the HDF5 library.

    <b>Build option</b>

    Support for HDF5 is optional. It is included when the \c BUILD_HDF5 build option is enabled,
    which requires the HDF5 C library to be installed on the system. Otherwise the code is built
    without any dependency on HDF5. In that case, the available() function returns false, and the
    functions that open a suite or a bundle throw a fatal error.

    <b>Concepts</b>

    An HDF5 file is organized like a small file system: it has a root group, which can contain
    further groups as well as datasets, the actual blocks of stored data. Both groups and datasets
    can carry attributes, small named pieces of metadata. The wrapper represents the following
    concepts, each by its own class:

    - A \em bundle (H5BundleR) is the counterpart of a regular input file. It is an HDF5 group
      containing one or more datasets, with names fixed by the SKIRT data model for the type of
      file. For example, a bundle that corresponds to a text column file contains a dataset for
      each column.

    - A \em dataset (H5DatasetR) is a rectangular, multidimensional array of values of the same
      type, located in a bundle.

    - A \em suite (H5SuiteR) is an HDF5 group that holds a set of bundles, for example the input
      files for a single simulation. This allows storing the input for several simulations in a
      single HDF5 file. A suite can be the root group of the file, or a group at any nesting level.

    The names of these classes are suffixed with "R" (for read-only), so that future extensions
    can differentiate between reading and writing.

    <b>Paths</b>

    None of the classes represents an open file; managing the open files is an internal concern
    of the wrapper. Instead, the functions of this class take a path that identifies both the HDF5
    file and the suite or bundle inside it, in the form

        <filepath>:<suite>/<bundle>

    where <tt>\<filepath\></tt> is the path of the HDF5 file, which must have the <tt>.hdf5</tt>
    filename extension, and the portion after the colon is the path of the group inside the file,
    with nested groups separated by forward slashes. The suite portion can be empty, in which case
    the bundle is located in the root group of the file, or it can consist of several nested
    groups. For example:

        /data/campaign7.hdf5:galaxy042/stars.txt

    The colon that separates the two portions is the one immediately following the first
    occurrence of the <tt>.hdf5</tt> filename extension.

    Combining the name of an input file with the input path and the HDF5 file specified by the
    user to construct such a path, and deciding whether to use a regular file or a bundle, is
    outside the scope of the wrapper. This is handled by the FilePaths class and its clients.

    <b>Object lifetime and thread safety</b>

    Each object handed out by the wrapper acquires its underlying HDF5 resource on construction
    and releases it in its destructor. The objects can be moved but not copied. An object remains
    fully valid after the object from which it was obtained has been destroyed. See the H5ObjectR
    class for more information.

    The HDF5 library cannot be assumed to be thread-safe on its own, since that depends on how
    the installed library was built. All calls into the library are therefore protected by a
    single mutex internal to the wrapper, so that all functions of the wrapper are thread-safe.
    Any number of independent processes may open and read the same HDF5 file at the same time. */
class H5Lib final
{
public:
    /** This function returns true if the code was built with HDF5 support, and false otherwise. */
    static bool available();

    /** This function opens the suite at the specified path for reading. The path has the form
        <tt>\<filepath\>:\<suite\></tt> as described in the class header; if the suite portion is
        empty (with or without the colon), the function opens the root group of the file. The
        function throws a fatal error if the code was built without HDF5 support, if the path is
        malformed, or if the file or the suite cannot be opened. */
    static H5SuiteR openSuite(const string& path);

    /** This function opens the bundle at the specified path for reading, without going through
        its suite first. The path has the form <tt>\<filepath\>:\<suite\>/\<bundle\></tt> as
        described in the class header. The function throws a fatal error if the code was built
        without HDF5 support, if the path is malformed, or if the file or the bundle cannot be
        opened. */
    static H5BundleR openBundle(const string& path);

    /** This function returns true if there is a bundle at the specified path, i.e. if the HDF5
        file exists and holds a group at the location of the bundle, and false if the file does
        not exist or does not hold such a group. The path has the same form as for the
        openBundle() function. The function throws a fatal error if the code was built without
        HDF5 support, if the path is malformed, or if the file exists but cannot be opened as an
        HDF5 file. It allows a caller to report a missing input file in its own words, for example
        to mention a regular file that is missing as well. */
    static bool hasBundle(const string& path);
};

////////////////////////////////////////////////////////////////////

#endif
